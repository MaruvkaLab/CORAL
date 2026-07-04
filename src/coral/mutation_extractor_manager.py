import os
import gzip
import json
import re
import multiprocessing
import subprocess
from io import TextIOWrapper
from collections import defaultdict

try:
    from .utils import log
except ImportError:
    def log(message, verbose=True):
        if verbose:
            print(message)

import pandas as pd

CHR_IDX, POSITION_IDX, REF_NUC_IDX, N_READS_1_IDX, NUC_1_IDX, N_READS_2_IDX, NUC_2_IDX = range(7)
PREV_IDX, CUR_IDX, NEXT_IDX = 0, 1, 2
REF_IDX, TAXA1_IDX, TAXA2_IDX = 0, 1, 2
FLANK = 2
MIN_DEPTH = 1  # Minimum depth threshold for quality check

# Matches read-start (^ + mapQ char), read-end ($), or an indel marker (+N / -N).
# Digits are captured so the N following indel bases can be skipped in one pass.
_CLEAN_RE = re.compile(r'\^.|\$|[+-](\d*)')


def clean_bases(s):
    """Strip read-start/end markers and indel notation from an mpileup bases field.

    Fast path: every marker this removes (^x, $, +N…, -N…) contains one of the
    four characters ^ $ + - . When none are present — the overwhelmingly common
    case for CORAL's shallow pileups — the regex below matches nothing and would
    return s unchanged, so we skip it entirely. This is byte-for-byte identical
    to the regex path but avoids the regex-engine call and list building.
    """
    if '^' not in s and '$' not in s and '+' not in s and '-' not in s:
        return s
    out = []
    pos = 0
    for m in _CLEAN_RE.finditer(s):
        if m.start() < pos:
            continue  # falls inside a just-skipped indel run; can't happen in valid pileup syntax
        out.append(s[pos:m.start()])
        digits = m.group(1)
        if digits is not None:
            pos = m.end() + (int(digits) if digits else 0)
        else:
            pos = m.end()
    out.append(s[pos:])
    return ''.join(out)


class PileupLine:
    """Wraps one parsed pileup row. clean_bases() has a cheap fast path (see
    above), so re-cleaning a field costs a few membership tests; we skip the
    per-line cache dict that this used to allocate on every one of ~10^8 lines."""
    __slots__ = ('fields',)

    def __init__(self, fields):
        self.fields = fields

    def __getitem__(self, idx):
        return self.fields[idx]

    def clean(self, idx):
        return clean_bases(self.fields[idx])

    def nuc(self, idx):
        cleaned = clean_bases(self.fields[idx])
        return cleaned[0].upper() if cleaned else 'N'


def parse_line(line):
    parts = line.strip().split('\t')
    if len(parts) < 9:
        return None
    return PileupLine(parts[:5] + parts[6:-1])


def all_same(seq):
    return len(seq) > 0 and all(ch == seq[0] for ch in seq)


def quality_check(line):
    if line is None:
        return False
    fields = line.fields
    if '*' in fields[NUC_1_IDX] or '*' in fields[NUC_2_IDX]: # deletions
        return False
    if '+' in fields[NUC_1_IDX] or '+' in fields[NUC_2_IDX]: # insertions
        return False
    if int(fields[N_READS_1_IDX]) < MIN_DEPTH or int(fields[N_READS_2_IDX]) < MIN_DEPTH:
        return False
    nuc1 = line.clean(NUC_1_IDX).replace(',', '.').lower()
    nuc2 = line.clean(NUC_2_IDX).replace(',', '.').lower()
    return all_same(nuc1) and all_same(nuc2)


def consecutive(*lines):
    chrom = lines[0][CHR_IDX]
    positions = [int(line[POSITION_IDX]) for line in lines]
    return (all(line[CHR_IDX] == chrom for line in lines)
            and all(positions[i] + 1 == positions[i + 1] for i in range(len(positions) - 1)))


def extract_context(lines):
    sequences = [[], [], []]
    for line in lines:
        ref_nuc = line.nuc(REF_NUC_IDX)
        nuc1 = line.nuc(NUC_1_IDX)
        nuc2 = line.nuc(NUC_2_IDX)
        sequences[REF_IDX].append(ref_nuc)
        sequences[TAXA1_IDX].append(ref_nuc if nuc1 in {',', '.'} else nuc1)
        sequences[TAXA2_IDX].append(ref_nuc if nuc2 in {',', '.'} else nuc2)
    return sequences


def detect_mutation_triplet(triplets):
    t1_mut = t2_mut = t1_3mer = t2_3mer = None

    if triplets[REF_IDX][PREV_IDX] == triplets[TAXA1_IDX][PREV_IDX] == triplets[TAXA2_IDX][PREV_IDX] and \
       triplets[REF_IDX][NEXT_IDX] == triplets[TAXA1_IDX][NEXT_IDX] == triplets[TAXA2_IDX][NEXT_IDX]:

        ref_base = triplets[REF_IDX][CUR_IDX]
        context = triplets[REF_IDX][PREV_IDX] + ref_base + triplets[REF_IDX][NEXT_IDX]
        if ref_base == triplets[TAXA1_IDX][CUR_IDX] and ref_base != triplets[TAXA2_IDX][CUR_IDX]:
            t2_mut = f"{context[0]}[{ref_base}>{triplets[TAXA2_IDX][CUR_IDX]}]{context[2]}"
            t1_3mer, t2_3mer = context, context
        elif ref_base == triplets[TAXA2_IDX][CUR_IDX] and ref_base != triplets[TAXA1_IDX][CUR_IDX]:
            t1_mut = f"{context[0]}[{ref_base}>{triplets[TAXA1_IDX][CUR_IDX]}]{context[2]}"
            t1_3mer, t2_3mer = context, context
        elif ref_base == triplets[TAXA2_IDX][CUR_IDX] == triplets[TAXA1_IDX][CUR_IDX]:
            t1_3mer, t2_3mer = context, context

    return t1_mut, t2_mut, t1_3mer, t2_3mer


def scan_pileup(line_iter, on_mut1=None, on_mut2=None):
    """Single linear scan over pileup lines with a 3-position sliding window.

    Returns (mut1, mut2, triplet1, triplet2) count dicts. For each detected
    mutation calls on_mut1/on_mut2(chrom, pos, mutation) when provided — the
    serial path uses these to stream CSV rows, the parallel path to collect
    them. This is the single source of truth for the scan: both the serial
    MutationExtractor and the parallel driver call it, so their results cannot
    diverge.
    """
    mut1 = defaultdict(int)
    mut2 = defaultdict(int)
    trip1 = defaultdict(int)
    trip2 = defaultdict(int)

    it = iter(line_iter)
    first = next(it, None)
    second = next(it, None)
    if first is None or second is None:
        return mut1, mut2, trip1, trip2

    line_fields = [None, parse_line(first), parse_line(second)]
    qc_flags = [False, quality_check(line_fields[1]), quality_check(line_fields[2])]

    for line in it:
        line_fields = [line_fields[1], line_fields[2], parse_line(line)]
        qc_flags = [qc_flags[1], qc_flags[2], quality_check(line_fields[2])]

        if all(qc_flags) and consecutive(*line_fields):
            triplets = extract_context(line_fields)
            m1, m2, t1, t2 = detect_mutation_triplet(triplets)
            chrom = line_fields[1][CHR_IDX]
            pos = int(line_fields[1][POSITION_IDX])

            if m1:
                mut1[m1] += 1
                if on_mut1:
                    on_mut1(chrom, pos, m1)
            if m2:
                mut2[m2] += 1
                if on_mut2:
                    on_mut2(chrom, pos, m2)
            if t1:
                trip1[t1] += 1
            if t2:
                trip2[t2] += 1

    return mut1, mut2, trip1, trip2


def _write_extractor_outputs(mut1, mut2, trip1, trip2, out_json1, out_json2, trip_json1, trip_json2):
    with open(out_json1, 'w') as f:
        json.dump(mut1, f, indent=2)
    with open(out_json2, 'w') as f:
        json.dump(mut2, f, indent=2)
    with open(trip_json1, 'w') as f:
        json.dump(trip1, f, indent=2)
    with open(trip_json2, 'w') as f:
        json.dump(trip2, f, indent=2)


class MutationExtractor:
    def __init__(self, reference, taxon1, taxon2, pileup_file, mutation_output_dir, triplet_output_dir,
                 no_full_mutations=False, no_cache=False, verbose=True):
        self.reference = reference
        self.taxon1 = taxon1
        self.taxon2 = taxon2
        self.pileup_file = pileup_file
        self.mutation_output_dir = mutation_output_dir
        self.triplet_output_dir = triplet_output_dir
        self.no_full_mutations = no_full_mutations
        self.no_cache = no_cache
        self.verbose = verbose

        self.out_json1 = os.path.join(self.mutation_output_dir, f"{taxon1}__{taxon2}__{reference}__mutations.json")
        self.out_json2 = os.path.join(self.mutation_output_dir, f"{taxon2}__{taxon1}__{reference}__mutations.json")
        self.trip_out_json1 = os.path.join(self.triplet_output_dir, f"{self.taxon1}__{self.taxon2}__{self.reference}__triplets.json")
        self.trip_out_json2 = os.path.join(self.triplet_output_dir, f"{self.taxon2}__{self.taxon1}__{self.reference}__triplets.json")

        self.csv_path1 = None if no_full_mutations else os.path.join(self.mutation_output_dir, f"{taxon1}__{taxon2}__{reference}__mutations.csv.gz")
        self.csv_path2 = None if no_full_mutations else os.path.join(self.mutation_output_dir, f"{taxon2}__{taxon1}__{reference}__mutations.csv.gz")

    def extract(self):
        os.makedirs(self.mutation_output_dir, exist_ok=True)
        os.makedirs(self.triplet_output_dir, exist_ok=True)

        jsons_exist = all(os.path.exists(p) for p in [self.out_json1, self.out_json2, self.trip_out_json1, self.trip_out_json2])
        csvs_exist = (self.no_full_mutations or
                      all(os.path.exists(p) for p in [self.csv_path1, self.csv_path2]))

        if not self.no_cache and jsons_exist and csvs_exist:
            log("Mutation counts already exist. Skipping.", self.verbose)
            return

        csv1 = csv2 = None
        on_mut1 = on_mut2 = None
        if not self.no_full_mutations:
            header = "chromosome,position,mutation\n"
            csv1 = gzip.open(self.csv_path1, 'wt')
            csv1.write(header)
            csv2 = gzip.open(self.csv_path2, 'wt')
            csv2.write(header)
            on_mut1 = lambda chrom, pos, m: csv1.write(f"{chrom},{pos},{m}\n")
            on_mut2 = lambda chrom, pos, m: csv2.write(f"{chrom},{pos},{m}\n")

        with gzip.open(self.pileup_file, 'rt') as f:
            species_mut1, species_mut2, species_triplet1, species_triplet2 = scan_pileup(f, on_mut1, on_mut2)

        if csv1:
            csv1.close()
        if csv2:
            csv2.close()

        _write_extractor_outputs(species_mut1, species_mut2, species_triplet1, species_triplet2,
                                 self.out_json1, self.out_json2, self.trip_out_json1, self.trip_out_json2)

        log(f"Saved mutation counts to {self.out_json1} and {self.out_json2}", self.verbose)
        log(f"Saved triplet counts to {self.trip_out_json1} and {self.trip_out_json2}", self.verbose)

    def detect_mutation_triplet(self, triplets):
        return detect_mutation_triplet(triplets)


def _read_fai_chroms(fai_path):
    """Return [(chrom, length), ...] in reference (.fai) order."""
    chroms = []
    with open(fai_path) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 2:
                chroms.append((parts[0], int(parts[1])))
    return chroms


def _chroms_with_reads(bams):
    """Union of reference names with >=1 mapped read in any BAM, via
    `samtools idxstats` (reads only the index, near-instant). Chromosomes absent
    here produce no pileup lines, so skipping them cannot change any output -- it
    just avoids launching an mpileup for empty scaffolds."""
    have = set()
    for bam in bams:
        out = subprocess.run(["samtools", "idxstats", bam],
                             capture_output=True, text=True, check=True).stdout
        for line in out.splitlines():
            parts = line.split('\t')
            if len(parts) >= 3 and parts[0] != '*' and int(parts[2]) > 0:
                have.add(parts[0])
    return have


def _extract_region(chrom, ref_fasta, bams, no_full_mutations):
    """Worker: generate this chromosome's pileup via index-based
    `samtools mpileup -r <chrom>` (so no worker ever streams the whole file),
    then scan it with the shared scan_pileup. The mpileup options match
    Pileup.generate exactly (-B -d 100), so a chromosome's output is byte-
    identical to its slice of the whole-genome pileup. Returns (chrom, count
    dicts, CSV rows)."""
    rows1 = []
    rows2 = []
    on1 = on2 = None
    if not no_full_mutations:
        on1 = lambda c, pos, m: rows1.append(f"{c},{pos},{m}\n")
        on2 = lambda c, pos, m: rows2.append(f"{c},{pos},{m}\n")

    cmd = ["samtools", "mpileup", "-f", ref_fasta, "-B", "-d", "100", "-r", chrom] + list(bams)
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    try:
        with TextIOWrapper(proc.stdout) as stream:
            mut1, mut2, trip1, trip2 = scan_pileup(stream, on1, on2)
    finally:
        proc.wait()
    return chrom, dict(mut1), dict(mut2), dict(trip1), dict(trip2), rows1, rows2


class ParallelMutationExtractor:
    """Chromosome-parallel drop-in for MutationExtractor producing identical
    output files. One task per chromosome generates that chromosome's pileup
    directly from the indexed BAMs (`samtools mpileup -r`) and scans it, so no
    worker reads the whole genome; a pool runs them concurrently and results are
    merged (counts summed, CSV rows concatenated in reference/.fai order). Only
    chromosomes that actually carry reads are processed."""

    def __init__(self, reference, taxon1, taxon2, ref_fasta, bams, mutation_output_dir, triplet_output_dir,
                 fai_path, cores, no_full_mutations=False, no_cache=False, verbose=True):
        self.reference = reference
        self.taxon1 = taxon1
        self.taxon2 = taxon2
        self.ref_fasta = ref_fasta
        self.bams = list(bams)
        self.mutation_output_dir = mutation_output_dir
        self.triplet_output_dir = triplet_output_dir
        self.fai_path = fai_path
        self.cores = cores
        self.no_full_mutations = no_full_mutations
        self.no_cache = no_cache
        self.verbose = verbose

        self.out_json1 = os.path.join(mutation_output_dir, f"{taxon1}__{taxon2}__{reference}__mutations.json")
        self.out_json2 = os.path.join(mutation_output_dir, f"{taxon2}__{taxon1}__{reference}__mutations.json")
        self.trip_out_json1 = os.path.join(triplet_output_dir, f"{taxon1}__{taxon2}__{reference}__triplets.json")
        self.trip_out_json2 = os.path.join(triplet_output_dir, f"{taxon2}__{taxon1}__{reference}__triplets.json")
        self.csv_path1 = None if no_full_mutations else os.path.join(mutation_output_dir, f"{taxon1}__{taxon2}__{reference}__mutations.csv.gz")
        self.csv_path2 = None if no_full_mutations else os.path.join(mutation_output_dir, f"{taxon2}__{taxon1}__{reference}__mutations.csv.gz")

    def extract(self):
        os.makedirs(self.mutation_output_dir, exist_ok=True)
        os.makedirs(self.triplet_output_dir, exist_ok=True)

        jsons_exist = all(os.path.exists(p) for p in [self.out_json1, self.out_json2, self.trip_out_json1, self.trip_out_json2])
        csvs_exist = (self.no_full_mutations or all(os.path.exists(p) for p in [self.csv_path1, self.csv_path2]))
        if not self.no_cache and jsons_exist and csvs_exist:
            log("Mutation counts already exist. Skipping.", self.verbose)
            return

        chrom_lengths = _read_fai_chroms(self.fai_path)
        chrom_order = [c for c, _ in chrom_lengths]
        with_reads = _chroms_with_reads(self.bams)
        # one task per chromosome that carries reads, longest-first so the big
        # ones start early and dominate the makespan; empty scaffolds are skipped.
        tasks = [c for c, _ in sorted(chrom_lengths, key=lambda cl: cl[1], reverse=True) if c in with_reads]
        n_workers = max(1, min(self.cores, len(tasks) or 1))
        log(f"Extracting mutations in parallel: {len(tasks)} chromosomes with reads over {n_workers} workers...", self.verbose)

        args = [(c, self.ref_fasta, self.bams, self.no_full_mutations) for c in tasks]
        if not args:
            results = []
        elif n_workers <= 1:
            results = [_extract_region(*a) for a in args]
        else:
            # chunksize=1: dispatch one chromosome per grab so the few big arms
            # land on different workers. The default chunksize bundles contiguous
            # (here longest-first) tasks, which would pile every big chromosome
            # onto a single worker and serialize the heavy work.
            with multiprocessing.Pool(n_workers) as pool:
                results = pool.starmap(_extract_region, args, chunksize=1)

        mut1 = defaultdict(int)
        mut2 = defaultdict(int)
        trip1 = defaultdict(int)
        trip2 = defaultdict(int)
        rows1 = {}
        rows2 = {}
        for chrom, m1, m2, t1, t2, r1, r2 in results:
            for k, v in m1.items():
                mut1[k] += v
            for k, v in m2.items():
                mut2[k] += v
            for k, v in t1.items():
                trip1[k] += v
            for k, v in t2.items():
                trip2[k] += v
            if r1:
                rows1[chrom] = r1
            if r2:
                rows2[chrom] = r2

        if not self.no_full_mutations:
            header = "chromosome,position,mutation\n"
            with gzip.open(self.csv_path1, 'wt') as c1:
                c1.write(header)
                for chrom in chrom_order:
                    for row in rows1.get(chrom, ()):
                        c1.write(row)
            with gzip.open(self.csv_path2, 'wt') as c2:
                c2.write(header)
                for chrom in chrom_order:
                    for row in rows2.get(chrom, ()):
                        c2.write(row)

        _write_extractor_outputs(mut1, mut2, trip1, trip2,
                                 self.out_json1, self.out_json2, self.trip_out_json1, self.trip_out_json2)
        log(f"Saved mutation counts to {self.out_json1} and {self.out_json2}", self.verbose)
        log(f"Saved triplet counts to {self.trip_out_json1} and {self.trip_out_json2}", self.verbose)


class FiveMerExtractor:
    def __init__(self, reference, taxon1, taxon2, pileup_file, output_dir, no_cache=False, verbose=True):
        self.reference = reference
        self.taxon1 = taxon1
        self.taxon2 = taxon2
        self.pileup_file = pileup_file
        self.json1_path = os.path.join(output_dir, f"{taxon1}__{taxon2}__{reference}__5mers.json")
        self.json2_path = os.path.join(output_dir, f"{taxon2}__{taxon1}__{reference}__5mers.json")
        self.no_cache = no_cache
        self.verbose = verbose
        os.makedirs(output_dir, exist_ok=True)

    def detect_mutation_5mer(self, five_mers):
        flank_indices = [0, 1, 3, 4]
        center = 2
        for i in flank_indices:
            if not (five_mers[REF_IDX][i] == five_mers[TAXA1_IDX][i] == five_mers[TAXA2_IDX][i]):
                return None, None
        ref_base = five_mers[REF_IDX][center]
        t1_base = five_mers[TAXA1_IDX][center]
        t2_base = five_mers[TAXA2_IDX][center]
        context = ''.join(five_mers[REF_IDX])
        t1_mut = f"{context[:2]}[{ref_base}>{t1_base}]{context[3:]}" if t1_base != ref_base and t2_base == ref_base else None
        t2_mut = f"{context[:2]}[{ref_base}>{t2_base}]{context[3:]}" if t2_base != ref_base and t1_base == ref_base else None
        return t1_mut, t2_mut

    def extract(self):
        if all(os.path.exists(p) for p in [self.json1_path, self.json2_path]) and not self.no_cache:
            log("5-mer mutation files exist. Skipping.", self.verbose)
            return self.json1_path, self.json2_path

        species_mut1 = defaultdict(int)
        species_mut2 = defaultdict(int)

        with gzip.open(self.pileup_file, 'rt') as f:
            line_fields = [None] * (2 * FLANK + 1)
            qc_flags = [False] * (2 * FLANK + 1)
            for i in range(2 * FLANK):
                line_fields[i] = parse_line(f.readline())
                qc_flags[i] = quality_check(line_fields[i])

            for line in f:
                line_fields = line_fields[1:] + [parse_line(line)]
                qc_flags = qc_flags[1:] + [quality_check(line_fields[-1])]

                if all(qc_flags) and consecutive(*line_fields):
                    five_mers = extract_context(line_fields)
                    mut1, mut2 = self.detect_mutation_5mer(five_mers)
                    if mut1:
                        species_mut1[mut1] += 1
                    if mut2:
                        species_mut2[mut2] += 1

        with open(self.json1_path, 'w') as f:
            json.dump(species_mut1, f, indent=2)
        with open(self.json2_path, 'w') as f:
            json.dump(species_mut2, f, indent=2)

        log(f"Written: {self.json1_path}, {self.json2_path}", self.verbose)

        return self.json1_path, self.json2_path


class TripletExtractor:
    def __init__(self, reference, taxon1, taxon2, pileup_file, output_dir, no_cache=False, verbose=True):
        self.reference = reference
        self.taxon1 = taxon1
        self.taxon2 = taxon2
        self.pileup_file = pileup_file
        self.output_dir = output_dir
        self.no_cache = no_cache
        self.verbose = verbose

        os.makedirs(self.output_dir, exist_ok=True)

        self.out_json1 = os.path.join(
            self.output_dir, f"{self.taxon1}__{self.taxon2}__{self.reference}__triplets.json"
        )
        self.out_json2 = os.path.join(
            self.output_dir, f"{self.taxon2}__{self.taxon1}__{self.reference}__triplets.json"
        )

    def relevant_triplet(self, triplets):
        if triplets[REF_IDX][PREV_IDX] == triplets[TAXA1_IDX][PREV_IDX] == triplets[TAXA2_IDX][PREV_IDX] and \
           triplets[REF_IDX][NEXT_IDX] == triplets[TAXA1_IDX][NEXT_IDX] == triplets[TAXA2_IDX][NEXT_IDX]:
            if triplets[REF_IDX][CUR_IDX] == triplets[TAXA1_IDX][CUR_IDX] or \
               triplets[REF_IDX][CUR_IDX] == triplets[TAXA2_IDX][CUR_IDX]:
                return True
        return False

    def extract(self):
        if all(os.path.exists(p) for p in [self.out_json1, self.out_json2]) and not self.no_cache:
            log("Triplet counts already exist. Skipping.", self.verbose)
            return

        triplet_dict1 = defaultdict(int)
        triplet_dict2 = defaultdict(int)

        with gzip.open(self.pileup_file, 'rt') as f:
            line_fields = [None, parse_line(f.readline()), parse_line(f.readline())]
            qc_flags = [False, quality_check(line_fields[1]), quality_check(line_fields[2])]

            for line in f:
                line_fields = [line_fields[1], line_fields[2], parse_line(line)]
                qc_flags = [qc_flags[1], qc_flags[2], quality_check(line_fields[2])]

                if all(qc_flags) and consecutive(*line_fields):
                    triplets = extract_context(line_fields)
                    if self.relevant_triplet(triplets):
                        triplet = ''.join(triplets[REF_IDX])
                        triplet_dict1[triplet] += 1
                        triplet_dict2[triplet] += 1

        with open(self.out_json1, 'w') as f:
            json.dump(triplet_dict1, f, indent=2)
        with open(self.out_json2, 'w') as f:
            json.dump(triplet_dict2, f, indent=2)

        log(f"Triplet dictionaries written to:\n  • {self.out_json1}\n  • {self.out_json2}", self.verbose)


class MutationNormalizer:
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    valid_bases = {'A', 'C', 'G', 'T'}
    mutation_pattern = re.compile(r"^[ACGT]\[[ACGT]>[ACGT]\][ACGT]$")

    def __init__(self, input_dir, output_dir=None, divergence_time=None, verbose=True):
        self.input_dir = input_dir
        self.mutation_dir = os.path.join(input_dir, "Mutations")
        self.triplet_dir = os.path.join(input_dir, "Triplets")
        self.output_dir = output_dir or os.path.join(input_dir, "Tables")
        self.divergence_time = divergence_time
        self.verbose = verbose
        os.makedirs(self.output_dir, exist_ok=True)

        self.all_collapsed_mut = {}
        self.all_norm_mut = {}
        self.all_scaled_mut = {}
        self.all_triplets = {}

    def load_json(self, path):
        with open(path, 'r') as f:
            return json.load(f)

    def save_json(self, obj, path):
        with open(path, 'w') as f:
            json.dump(obj, f, indent=2)

    def collapse_triplets(self, triplet_dict):
        collapsed = defaultdict(int)
        for triplet, count in triplet_dict.items():
            if triplet[1] in {'G', 'A'}:
                rc = ''.join(self.complement.get(nuc, nuc) for nuc in reversed(triplet))
                collapsed[rc] += count
            else:
                collapsed[triplet] += count
        return collapsed

    def get_complement(self, mutation):
        comp = [self.complement[nuc] if nuc in self.complement else nuc for nuc in mutation]
        comp[0], comp[-1] = comp[-1], comp[0]
        return ''.join(comp)

    def collapse_mutations(self, mutation_dict):
        collapsed = defaultdict(int)
        for mutation, count in mutation_dict.items():
            if mutation[2] in {'A', 'G'}:
                collapsed[self.get_complement(mutation)] += int(count)
            else:
                collapsed[mutation] += int(count)
        return collapsed

    def filter_mutations_dict(self, d):
        return {k: v for k, v in d.items() if self.mutation_pattern.match(k)}

    def filter_triplets_dict(self, d):
        return {k: v for k, v in d.items() if 'N' not in k and all(x in self.valid_bases for x in k)}

    def normalize_by_triplets(self, mutations, triplets):
        return {
            k: v / triplets.get(f"{k[0]}{k[2]}{k[-1]}", 1) if triplets.get(f"{k[0]}{k[2]}{k[-1]}", 0) > 0 else 0
            for k, v in mutations.items()
        }

    def scale_counts(self, d, target_sum=10000):
        total = sum(d.values())
        if total == 0:
            return {k: 0 for k in d}
        return {k: round(v / total * target_sum) for k, v in d.items()}

    def normalize(self):
        for file in os.listdir(self.mutation_dir):
            if not file.endswith(".json"):
                continue

            mutation_path = os.path.join(self.mutation_dir, file)
            triplet_path = os.path.join(self.triplet_dir, file.replace("mutations", "triplets"))
            key = file.replace('.json', '')

            if not os.path.exists(triplet_path):
                continue

            mutations = self.filter_mutations_dict(self.load_json(mutation_path))
            triplets = self.filter_triplets_dict(self.load_json(triplet_path))

            collapsed_mut = self.collapse_mutations(mutations)
            collapsed_tri = self.collapse_triplets(triplets)

            self.all_collapsed_mut[key] = collapsed_mut
            self.all_triplets[key] = collapsed_tri

            norm_mut = self.normalize_by_triplets(collapsed_mut, collapsed_tri)
            scaled_mut = self.scale_counts(collapsed_mut)
            scaled_norm_mut = self.scale_counts(norm_mut)

            self.all_norm_mut[key] = scaled_norm_mut
            self.all_scaled_mut[key] = scaled_mut

        self._export()

    def _export(self):
        collapsed_mutations_df = pd.DataFrame(self.all_collapsed_mut)
        collapsed_mutations_df.to_csv(os.path.join(self.output_dir, "collapsed_mutations.tsv"), sep='\t')
        pd.DataFrame(self.all_norm_mut).to_csv(os.path.join(self.output_dir, "normalized_scaled.tsv"), sep='\t')
        pd.DataFrame(self.all_scaled_mut).to_csv(os.path.join(self.output_dir, "scaled_raw.tsv"), sep='\t')
        triplets_df = pd.DataFrame(self.all_triplets)
        triplets_df.to_csv(os.path.join(self.output_dir, "triplets.tsv"), sep='\t')

        mutations_per_triplet = collapsed_mutations_df.sum() / triplets_df.sum()
        log("Mutations per triplet from divergence:", self.verbose)
        for col in mutations_per_triplet.index:
            log(f"{col.split('__')[0]}: {mutations_per_triplet[col]:.2e}", self.verbose)

        if self.divergence_time:
            mutation_rates = mutations_per_triplet / (float(self.divergence_time) * 1_000_000)
            log("Estimated mutation rates per site per year:", self.verbose)
            for col in mutation_rates.index:
                log(f"{col.split('__')[0]}: {mutation_rates[col]:.2e}", self.verbose)
