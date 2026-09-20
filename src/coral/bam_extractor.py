"""Pair-mode calls made straight from the final BAM, with no pileup in between.

The pileup path spends its time on one text line per covered base; the reads are
~75x fewer. Here each tile of the reference is walked read by read into a handful
of per-position arrays, and NumPy finds the clean 3-base windows and the calls.
Alignment is processed in tiles, so memory stays flat on any genome.

This reproduces `samtools mpileup -B -d 0` piped into `scan_pair`, which stays in
`mutation_extractor_manager` as the reference implementation:

  * reads with UNMAP / SECONDARY / QCFAIL / DUP are skipped
  * a deletion flags the deleted positions (mpileup writes `*` there), an
    insertion flags the base *before* it (`+N`), and `-N` is not an indel event;

`scan_region` and `tile_spans`/`run_tiles` are deliberately mode-independent: the
CIGAR walk and the tiling are the same job for any number of samples.
"""

import gzip
import json
import multiprocessing
import os
from collections import defaultdict

import numpy as np
import pysam

from .mutation_extractor_manager import _chroms_with_reads, _read_fai_chroms, fold
from .utils import log

# mpileup's default --ff: UNMAP, SECONDARY, QCFAIL, DUP.
SKIP_FLAGS = 0x4 | 0x100 | 0x200 | 0x400

CIGAR_ALIGNED = (0, 7, 8)   # M, =, X
CIGAR_INSERT = 1
CIGAR_DELETE = 2
CIGAR_REFSKIP = 3
CIGAR_SOFTCLIP = 4

_B = "ACGT"
_CODE = np.full(256, 4, np.uint8)       # ASCII -> 0..3, 4 for anything else (N, ambiguity codes)
for _i, _b in enumerate(_B):
    _CODE[ord(_b)] = _i

# A 3-mer is packed as (b0 << 4) | (b1 << 2) | b2, so a (reference, target) pair
# of 3-mers indexes a flat 64x64 table. Every called window is one lookup.
_TRIPLET = [_B[c >> 4] + _B[(c >> 2) & 3] + _B[c & 3] for c in range(64)]
_PAIR_IDX = np.full(64 * 64, -1, np.int32)
_CHANGE, _MUTATION = [], []
for _r in range(64):
    for _t in range(64):
        _ref3, _tgt3 = _TRIPLET[_r], _TRIPLET[_t]
        if _ref3[0] == _tgt3[0] and _ref3[2] == _tgt3[2] and _ref3[1] != _tgt3[1]:
            _PAIR_IDX[_r * 64 + _t] = len(_CHANGE)
            _CHANGE.append(f"{_ref3[0]}[{_ref3[1]}>{_tgt3[1]}]{_ref3[2]}")
            _MUTATION.append(fold(_ref3, _tgt3))


_OPEN = {}


def _opened(kind, path):
    """Per-process handle cache: a worker scans hundreds of tiles, and reopening
    the BAM index or the .fai for each one costs more than the scan."""
    handle = _OPEN.get((kind, path))
    if handle is None:
        handle = pysam.FastaFile(path) if kind == 'fasta' else pysam.AlignmentFile(path)
        _OPEN[(kind, path)] = handle
    return handle


def scan_region(bam, chrom, lo, hi):
    """Walk every read overlapping the 0-based half-open span [lo, hi) into five
    per-position arrays, one slot per reference base:

      covered   a read covers this position (mpileup would emit a line: depth >= 1)
      base      the first read base seen here, 0 where none
      disagree  two reads showed different bases
      deleted   a read's deletion covers this position -- mpileup's '*'
      inserted  a read has an insertion right after this position -- mpileup's '+N'
    """
    n = max(0, hi - lo)
    covered = np.zeros(n, bool)
    base = np.zeros(n, np.uint8)
    disagree = np.zeros(n, bool)
    deleted = np.zeros(n, bool)
    inserted = np.zeros(n, bool)
    if n == 0:
        return covered, base, disagree, deleted, inserted

    for read in bam.fetch(chrom, lo, hi):
        if read.flag & SKIP_FLAGS or read.cigartuples is None or read.query_sequence is None:
            continue
        seq = np.frombuffer(read.query_sequence.encode(), np.uint8)
        start = read.reference_start
        r = start
        q = 0
        for op, oplen in read.cigartuples:
            if op in CIGAR_ALIGNED:
                a, b = max(r, lo), min(r + oplen, hi)
                if a < b:
                    sl = slice(a - lo, b - lo)
                    s = seq[q + (a - r):q + (b - r)]
                    have = base[sl]
                    disagree[sl] |= (have != 0) & (have != s)
                    base[sl] = np.where(have == 0, s, have)
                    covered[sl] = True
                r += oplen
                q += oplen
            elif op == CIGAR_DELETE:
                a, b = max(r, lo), min(r + oplen, hi)
                if a < b:
                    sl = slice(a - lo, b - lo)
                    deleted[sl] = True
                    covered[sl] = True
                r += oplen
            elif op == CIGAR_INSERT:
                # '+N' hangs off the aligned base before the insertion, so an
                # insertion at the very start of a read has nowhere to land.
                if q > 0 and r > start and lo < r <= hi:
                    inserted[r - 1 - lo] = True
                q += oplen
            elif op == CIGAR_SOFTCLIP:
                q += oplen
            elif op == CIGAR_REFSKIP:
                r += oplen
    return covered, base, disagree, deleted, inserted


def tile_spans(fai_path, bam_path, tile_bp):
    """returns a list of tiles in reference (.fai) order, to be processed in chunks to avoid exploding memory."""
    lengths = _read_fai_chroms(fai_path)
    have = _chroms_with_reads([bam_path])
    return [(chrom, start, min(start + tile_bp, length), length)
            for chrom, length in lengths if chrom in have
            for start in range(0, length, tile_bp)]


def tile_scan_span(lo, hi, length, indel_window):
    """Computes the 0-based half-open span a tile must scan: its own [lo, hi) plus enough 
    padding for every middle base's flanks and indel neighbourhood, clipped to the contig."""
    pad = max(1, indel_window)
    return max(0, lo - pad), min(length, hi + pad)


def run_tiles(tasks, worker, cores, handle):
    """Run tasks (tiles) in parallel and pass each result to `handle` in task order."""
    if not tasks:
        return
    if not cores or cores < 2:
        for task in tasks:
            handle(worker(task))
        return
    batch = 4 * cores
    with multiprocessing.Pool(min(cores, len(tasks))) as pool:
        for i in range(0, len(tasks), batch):
            for result in pool.map(worker, tasks[i:i + batch], chunksize=1):
                handle(result)


class GapCounter:
    """Totals `non_consecutive` over tiles fed in reference order.

    Each tile counts the gap triples that lie wholly inside it and hands over its
    first and last two lines. A triple straddling a boundary is counted here, once,
    against the carried tail of the tiles before it -- so tiling cannot change the
    total, and a tile with no lines at all just passes the carry through.
    """

    def __init__(self):
        self.total = 0
        self._carry = []
        self._chrom = None

    def add(self, chrom, gaps, edge, tail, n_lines):
        if chrom != self._chrom:        # a triple never spans two chromosomes
            self._carry, self._chrom = [], chrom
        self.total += gaps
        seq = self._carry + edge
        for i in range(min(len(self._carry), max(0, len(seq) - 2))):
            (p0, c0), (_, c1), (p2, c2) = seq[i], seq[i + 1], seq[i + 2]
            if c0 and c1 and c2 and p2 != p0 + 2:
                self.total += 1
        self._carry = (self._carry + edge)[-2:] if n_lines <= 2 else tail


def _scan_pair_tile(task):
    """Processes one tile: the calls whose middle base lies in [lo, hi). It spans
    [lo - pad, hi + pad) so those middles have their flanks, and their indel
    neighbourhood, inside the tile; every middle belongs to exactly one tile."""
    bam_path, fasta_path, chrom, lo, hi, length, ref_mask, indel_window = task
    a, b = tile_scan_span(lo, hi, length, indel_window)
    n = b - a

    ref = np.frombuffer(_opened('fasta', fasta_path).fetch(chrom, a, b).upper().encode(), np.uint8)
    covered, base, disagree, deleted, inserted = scan_region(_opened('bam', bam_path), chrom, a, b)
    masked = ref_mask.mask_array(chrom, a + 1, b) if ref_mask is not None else None
    # The clean positions are the good positions, that passed QC to use for calling.
    clean = covered & ~deleted & ~inserted & ~disagree
    if masked is not None:
        clean &= ~masked
    
    inside = slice(lo - a, hi - a)
    cov = covered[inside]
    if masked is None:
        usable, lines_masked = cov, 0
    else:
        msk = masked[inside]
        usable, lines_masked = cov & ~msk, int((cov & msk).sum())
    lines_after_mask = int(usable.sum())

     # Counts why each unusable line was unusable.
    dl, ins, dis = deleted[inside], inserted[inside], disagree[inside]
    dropped = {'deletion': int((usable & dl).sum()),
               'insertion': int((usable & ~dl & ins).sum()),
               'reads_disagree': int((usable & ~dl & ~ins & dis).sum())}

    if n >= 3:
        mid = np.flatnonzero(clean[:-2] & clean[1:-1] & clean[2:]) + 1
        mid = mid[(mid >= lo - a) & (mid < hi - a)]
    else:
        mid = np.empty(0, np.int64)

    # scan_windows needs three adjacent *lines* clean, not three bases. A gap is
    # such a triple whose outer lines are not two apart.
    line_idx = np.flatnonzero(usable)
    line_clean = clean[inside][line_idx]
    if line_idx.size >= 3:
        triple = line_clean[:-2] & line_clean[1:-1] & line_clean[2:]
        consecutive = line_idx[2:] == line_idx[:-2] + 2
        gaps = int((triple & ~consecutive).sum())
    else:
        gaps = 0
    edge = [(int(lo + p), bool(c)) for p, c in
            zip(line_idx[:2].tolist(), line_clean[:2].tolist())]
    tail = [(int(lo + p), bool(c)) for p, c in
            zip(line_idx[-2:].tolist(), line_clean[-2:].tolist())]

    # Drop calls sitting within 'indel_window' bases of an indel. 
    n_near = 0
    if indel_window > 1 and mid.size:
        cs = np.concatenate(([0], np.cumsum((deleted | inserted).astype(np.int64))))
        left = np.maximum(mid - indel_window, 0)
        right = np.minimum(mid + indel_window + 1, n)
        near = (cs[right] - cs[left]) > 0
        n_near = int(near.sum())
        mid = mid[~near]

    rc, tc = _CODE[ref], _CODE[base]
    # Checks all bases are A/C/G/T, the flanks are identical, and whether the middle differs (mutation).
    acgt = ((rc[mid - 1] < 4) & (rc[mid] < 4) & (rc[mid + 1] < 4)
            & (tc[mid - 1] < 4) & (tc[mid] < 4) & (tc[mid + 1] < 4))
    flanks = (rc[mid - 1] == tc[mid - 1]) & (rc[mid + 1] == tc[mid + 1])
    same = rc[mid] == tc[mid]
    called = acgt & flanks & ~same
    kept = acgt & flanks
    classes = {'not_acgt': int((~acgt).sum()),
               'flanks_not_conserved': int((acgt & ~flanks).sum()),
               'identical': int((kept & same).sum()),
               'differs': int(called.sum())}

    def code3(codes, pos):
        return (codes[pos - 1].astype(np.int64) << 4) | (codes[pos].astype(np.int64) << 2) | codes[pos + 1]

    trips = {}
    if kept.any():
        ref_codes = code3(rc, mid[kept])
        for code, count in enumerate(np.bincount(ref_codes, minlength=64).tolist()):
            if count:
                trips[_TRIPLET[code]] = count

    muts, rows = {}, ''
    if called.any():
        d = mid[called]
        pair = _PAIR_IDX[code3(rc, d) * 64 + code3(tc, d)]
        for idx, count in enumerate(np.bincount(pair, minlength=len(_CHANGE)).tolist()):
            if count:
                muts[_MUTATION[idx]] = muts.get(_MUTATION[idx], 0) + count
        rows = ''.join(f"{chrom},{p + 1 + a},{_CHANGE[i]},{_MUTATION[i]}\n"
                       for p, i in zip(d.tolist(), pair.tolist()))

    return (chrom, rows, muts, trips, classes, dropped, lines_after_mask, lines_masked,
            gaps, edge, tail, int(line_idx.size), n_near)


class BamPairExtractor:
    """Pair mode: one target against the reference, scanned from the BAM.

    Writes exactly the files `PairExtractor` writes -- same names, the same CSV
    columns, the same mutation and triplet JSON, the same `run_summary.json` sections.
    """

    def __init__(self, reference, target, ref_fasta, bam, mutation_output_dir, triplet_output_dir,
                 cores=None, no_cache=False, verbose=True, ref_mask=None, indel_window=1,
                 tile_bp=1_000_000):
        self.reference = reference
        self.target = target
        self.ref_fasta = ref_fasta
        self.bam = bam
        self.cores = cores
        self.no_cache = no_cache
        self.verbose = verbose
        self.ref_mask = ref_mask
        self.indel_window = indel_window
        self.tile_bp = tile_bp
        self.mutation_output_dir = mutation_output_dir
        self.triplet_output_dir = triplet_output_dir

        self.out_json = os.path.join(mutation_output_dir, f"{target}__{reference}__mutations.json")
        self.csv_path = os.path.join(mutation_output_dir, f"{target}__{reference}__mutations.csv.gz")
        self.trip_out_json = os.path.join(triplet_output_dir, f"{target}__{reference}__triplets.json")
        self.classes_json = os.path.join(os.path.dirname(mutation_output_dir), "run_summary.json")

    def extract(self):
        os.makedirs(self.mutation_output_dir, exist_ok=True)
        os.makedirs(self.triplet_output_dir, exist_ok=True)
        if not self.no_cache and all(os.path.exists(p) for p in
                                     [self.out_json, self.csv_path, self.trip_out_json, self.classes_json]):
            log("Mutation counts already exist. Skipping.", self.verbose)
            return

        spans = tile_spans(self.ref_fasta + ".fai", self.bam, self.tile_bp)

        def tile_mask(chrom, lo, hi, length):
            """Just this tile's slice of the mask; None when it holds no repeats,
            so the worker skips mask_array and the masking entirely."""
            if not self.ref_mask:
                return None
            a, b = tile_scan_span(lo, hi, length, self.indel_window)
            sub = self.ref_mask.for_span(chrom, a + 1, b)   # the mask is 1-based inclusive
            return sub if sub else None

        tasks = [(self.bam, self.ref_fasta, chrom, lo, hi, length,
                  tile_mask(chrom, lo, hi, length), self.indel_window)
                 for chrom, lo, hi, length in spans]
        log(f"Scanning the BAM directly: {len(tasks)} tiles of {self.tile_bp:,} bp...", self.verbose)

        muts = defaultdict(int)
        trips = defaultdict(int)
        classes = defaultdict(int)
        dropped = defaultdict(int)
        lines = defaultdict(int)
        n_near = 0
        gap_counter = GapCounter()

        with gzip.open(self.csv_path, 'wt') as csv:
            csv.write("chromosome,position,change,mutation\n")

            def handle(result):
                nonlocal n_near
                (chrom, rows, m, t, cl, dr, after_mask, masked_lines,
                 gaps, edge, tail, n_lines, near) = result
                csv.write(rows)
                for dst, src in ((muts, m), (trips, t), (classes, cl), (dropped, dr)):
                    for k, v in src.items():
                        dst[k] += v
                lines['lines_after_mask'] += after_mask
                lines['lines_masked'] += masked_lines
                n_near += near
                gap_counter.add(chrom, gaps, edge, tail, n_lines)

            run_tiles(tasks, _scan_pair_tile, self.cores, handle)

        windows = {'non_consecutive': gap_counter.total}
        if self.indel_window > 1:
            windows['near_indel'] = n_near
        summary = {
            'pileup_lines': {k: v for k, v in dropped.items() if v} | dict(lines),
            'candidate_windows': windows,
            'site_classes': {k: v for k, v in classes.items() if v},
        }

        with open(self.out_json, 'w') as f:
            json.dump(dict(muts), f, indent=2)
        with open(self.trip_out_json, 'w') as f:
            json.dump(dict(trips), f, indent=2)
        with open(self.classes_json, 'w') as f:
            json.dump({section: dict(counts) for section, counts in summary.items()}, f, indent=2)
        log(f"Saved mutation counts to {self.out_json} and triplet counts to {self.trip_out_json}", self.verbose)
