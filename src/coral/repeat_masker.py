"""Repeat masking utilities for the CORAL pipeline.

This module runs an external repeat masker on a genome and stores the resulting
repeat intervals in a `RepeatMask`.

The mask is used to:
  * skip pseudo-reads whose fragment is mostly repetitive
    (`Genome.generate_fragment_fastq`);
  * skip pileup positions that fall inside repetitive regions of the outgroup.

Two repeat-masking backends are supported:
  * **windowmasker** (default) -- NCBI WindowMasker; de novo and
    k-mer-frequency based, with no species-specific repeat library required.
  * **repeatmasker** -- RepeatMasker; library-based and optionally configured
    for a specific species or clade.

Both backends are converted to the same interval representation, so downstream
code is independent of the masking tool used.
"""

import os
from array import array

import numpy as np

from .utils import run_cmd, log


def _merge(a):
    """Merge an (N, 2) array of 1-based inclusive intervals."""
    if len(a) == 0:
        return np.empty(0, np.int64), np.empty(0, np.int64)

    a = a[np.argsort(a[:, 0], kind="stable")]
    starts, ends = a.T

    # Maximum end reached up to each interval
    max_end = np.maximum.accumulate(ends)

    # Indices where a new merged interval begins
    new_group = np.r_[True, starts[1:] > max_end[:-1] + 1]
    idx = np.flatnonzero(new_group)

    return starts[idx], np.maximum.reduceat(ends, idx)


class RepeatMask:
    """Per-contig repeat intervals (1-based, inclusive) with fast lookup.

    Coordinates are 1-based inclusive to match RepeatMasker `.out`, samtools
    pileup positions, and CORAL's fragment read names -- so callers never convert.
    An empty mask (no file / masking disabled) is falsy and passes everything.
    """

    def __init__(self, intervals=None):
        # chrom -> (starts[], ends[]) parallel sorted arrays of merged intervals
        self._starts = {}
        self._ends = {}
        for chrom, ivs in (intervals or {}).items():
            s, e = _merge(np.asarray(ivs, dtype=np.int64).reshape(-1, 2))
            self._starts[chrom] = s
            self._ends[chrom] = e

    def __bool__(self):
        return any(len(v) for v in self._starts.values())

    @property
    def n_intervals(self):
        return sum(len(v) for v in self._starts.values())

    @property
    def n_masked_bases(self):
        return sum(int(self._ends[chrom].sum() - starts.sum()) + len(starts)
                   for chrom, starts in self._starts.items())

    def contains(self, chrom, pos):
        """Is 1-based position `pos` inside a repeat on `chrom`?"""
        starts = self._starts.get(chrom)
        if starts is None or len(starts) == 0:
            return False
        i = int(starts.searchsorted(pos, "right")) - 1
        return i >= 0 and pos <= self._ends[chrom][i]

    def for_contig(self, chrom):
        """A RepeatMask holding only `chrom`'s intervals (shares the underlying
        lists, no copy). Lets the parallel extractor ship each worker just its
        chromosome's mask instead of the whole genome's."""
        sub = RepeatMask()
        if chrom in self._starts:
            sub._starts[chrom] = self._starts[chrom]
            sub._ends[chrom] = self._ends[chrom]
        return sub

    def masked_bases(self, chrom, start, end):
        """Number of masked bases in the 1-based inclusive span [start, end]."""
        starts = self._starts.get(chrom)
        if starts is None or len(starts) == 0:
            return 0
        ends = self._ends[chrom]
        total = 0
        i = int(starts.searchsorted(end, "right")) - 1   # last interval starting <= end
        while i >= 0 and ends[i] >= start:
            total += min(end, int(ends[i])) - max(start, int(starts[i])) + 1
            i -= 1
        return total

    def masked_fraction(self, chrom, start, end):
        span = end - start + 1
        if span <= 0:
            return 0.0
        return self.masked_bases(chrom, start, end) / span

    # ---- construction ----

    @classmethod
    def from_out(cls, out_path):
        """Parse a RepeatMasker `.out` file into a RepeatMask.

        Columns (whitespace-delimited, 3 header lines): SW_score perc_div perc_del
        perc_ins QUERY QBEGIN QEND q_left strand repeat class/family ... -- so the
        query contig is field 5 and the 1-based inclusive span is fields 6-7."""
        intervals = {}
        with open(out_path) as fh:
            for line in fh:
                parts = line.split()
                if len(parts) < 7 or not parts[0].isdigit():
                    continue                          # header / blank
                chrom = parts[4]
                try:
                    start, end = int(parts[5]), int(parts[6])
                except ValueError:
                    continue
                if start > end:
                    start, end = end, start
                intervals.setdefault(chrom, array("q")).extend((start, end))
        return cls(intervals)

    @classmethod
    def from_windowmasker(cls, intervals_path):
        """Parse a WindowMasker `-outfmt interval` file into a RepeatMask.

        Format: a `>contig` header line, then `start - end` lines whose coordinates
        are 0-based inclusive -- so we shift both ends by +1 to our 1-based frame."""
        intervals = {}
        chrom = None
        with open(intervals_path) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    chrom = line[1:].split()[0]
                    continue
                a, _, b = line.partition('-')
                try:
                    start, end = int(a) + 1, int(b) + 1
                except ValueError:
                    continue
                if chrom is not None:
                    intervals.setdefault(chrom, array("q")).extend((start, end))
        return cls(intervals)

    @classmethod
    def from_bed(cls, bed_path):
        """Parse a BED (0-based, half-open) into a RepeatMask (1-based inclusive)."""
        intervals = {}
        with open(bed_path) as fh:
            for line in fh:
                if not line.strip() or line.startswith(('#', 'track', 'browser')):
                    continue
                c = line.split('\t')
                if len(c) < 3:
                    continue
                intervals.setdefault(c[0], array("q")).extend((int(c[1]) + 1, int(c[2])))
        return cls(intervals)

    def to_bed(self, bed_path):
        """Write merged intervals as BED (0-based, half-open) for inspection/cache."""
        with open(bed_path, 'w') as out:
            for chrom in sorted(self._starts):
                for s, e in zip(self._starts[chrom], self._ends[chrom]):
                    out.write(f"{chrom}\t{s - 1}\t{e}\n")


def run_repeatmasker(fasta_path, out_dir, species=None, cores=1,
                     no_cache=False, verbose=True):
    """Run RepeatMasker on `fasta_path`, writing outputs under `out_dir`; return the
    path to the `.out` annotation. Cached: skips the (slow) run if the `.out` already
    exists. `species` selects the RepeatMasker library clade (e.g. "drosophila")."""
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, os.path.basename(fasta_path) + ".out")

    if os.path.exists(out_path) and not no_cache:
        log(f"RepeatMasker output cached for {os.path.basename(fasta_path)}. Skipping.",
            verbose)
        return out_path

    # -pa is the number of parallel batch jobs; rmblast uses ~4 threads per job.
    pa = max(1, cores // 4)
    cmd = ["RepeatMasker", "-pa", str(pa), "-dir", out_dir]
    if species:
        cmd += ["-species", species]
    cmd += [fasta_path]

    log(f"Running RepeatMasker on {os.path.basename(fasta_path)} "
        f"(species={species or 'auto'}, pa={pa})...", verbose)
    run_cmd(cmd, verbose=verbose)

    if not os.path.exists(out_path):
        raise FileNotFoundError(
            f"RepeatMasker did not produce {out_path}; check its logs in {out_dir}")
    return out_path


def run_windowmasker(fasta_path, out_dir, dust=True, no_cache=False, verbose=True):
    """Run NCBI WindowMasker on `fasta_path` (two passes: build the k-mer count model,
    then emit masked intervals); return the path to the interval file. Cached: skips
    the run if the interval file already exists. `dust` also masks low-complexity."""
    os.makedirs(out_dir, exist_ok=True)
    base = os.path.join(out_dir, os.path.basename(fasta_path))
    counts_path = base + ".wm_counts"
    intervals_path = base + ".wm_intervals"

    if os.path.exists(intervals_path) and not no_cache:
        log(f"WindowMasker output cached for {os.path.basename(fasta_path)}. Skipping.",
            verbose)
        return intervals_path

    log(f"Running WindowMasker on {os.path.basename(fasta_path)}...", verbose)
    run_cmd(["windowmasker", "-in", fasta_path, "-infmt", "fasta",
             "-mk_counts", "-out", counts_path], verbose=verbose)
    run_cmd(["windowmasker", "-in", fasta_path, "-infmt", "fasta",
             "-ustat", counts_path, "-dust", "true" if dust else "false",
             "-outfmt", "interval", "-out", intervals_path], verbose=verbose)

    if not os.path.exists(intervals_path):
        raise FileNotFoundError(f"WindowMasker did not produce {intervals_path}")
    return intervals_path


def build_mask(fasta_path, out_dir, tool="windowmasker", species=None, cores=1,
               dust=True, no_cache=False, verbose=True):
    """Run the chosen masker (cached) and return a RepeatMask. Writes a merged BED
    alongside the raw output for inspection. `tool` is 'windowmasker' (default,
    library-free) or 'repeatmasker' (`species` selects its library clade)."""
    tool = (tool or "windowmasker").lower()
    if tool == "windowmasker":
        raw = run_windowmasker(fasta_path, out_dir, dust=dust,
                               no_cache=no_cache, verbose=verbose)
        mask = RepeatMask.from_windowmasker(raw)
    elif tool == "repeatmasker":
        raw = run_repeatmasker(fasta_path, out_dir, species=species, cores=cores,
                               no_cache=no_cache, verbose=verbose)
        mask = RepeatMask.from_out(raw)
    else:
        raise ValueError(f"unknown repeat masker tool {tool!r} "
                         "(expected 'windowmasker' or 'repeatmasker')")
    mask.to_bed(raw + ".repeats.bed")
    total_bp = _total_bp(fasta_path)
    masked = mask.n_masked_bases
    pct = f"{100 * masked / total_bp:.1f}% of {total_bp:,} bp" if total_bp else "genome length unknown"
    log(f"  repeat intervals: {mask.n_intervals}, masked bases: {masked:,} ({pct}, {tool})",
        verbose)
    return mask


def _total_bp(fasta_path):
    """Total sequence length of `fasta_path`, from its .fai if one exists."""
    fai = fasta_path + ".fai"
    if os.path.exists(fai):
        with open(fai) as fh:
            return sum(int(line.split("\t")[1]) for line in fh if line.strip())
    total = 0
    with open(fasta_path, "rb") as fh:
        for line in fh:
            if not line.startswith(b">"):
                total += len(line.strip())
    return total
