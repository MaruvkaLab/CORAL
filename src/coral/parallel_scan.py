"""Chromosome-parallel pileup scan for the multi-species extractor.

``MultipleSpeciesMutationExtractor.extract`` has two scan paths: the original
serial one, which streams a whole-genome ``samtools mpileup`` through a single
Python loop, and this one. Every bit of the detection logic is shared --
``_parse_line``, ``_quality_check``, ``_detect_site`` and ``_consecutive`` are
called on the extractor itself, never re-implemented here. What changes is only
*how the pileup is produced and traversed*: one task per chromosome generates
that chromosome's pileup directly from the indexed BAMs (``samtools mpileup
-r``) and scans it, so no worker ever reads the whole genome.

Why the split is exactly equivalent to the serial scan
------------------------------------------------------
1. ``_consecutive`` rejects any 3-line window that is not three consecutive
   positions on one chromosome, so a window straddling a chromosome boundary is
   already discarded by the serial scan. Splitting per chromosome can therefore
   neither create nor destroy a window.
2. The first two lines of a chromosome: in the serial scan, line 1 can only be a
   window centre together with a line from the *previous* chromosome (rejected),
   and line 2 gets exactly the window ``[1, 2, 3]`` -- which is precisely the
   first window the per-chromosome scan forms after seeding its buffer.
3. ``_quality_check`` and ``_detect_site`` look only at the window; there is no
   carried state.
4. mpileup emits ``3 + 3*(S-1)`` columns (the outgroup is ``-f``, not a BAM),
   which is what ``_parse_line``'s ``len(parts) >= n_species * 3`` guard expects;
   the BAM order stays ``species_list`` order, and that order is what defines the
   ``taxaK`` columns.
5. Chromosomes with no mapped reads produce no pileup lines at all, so skipping
   them cannot change the output.
6. mpileup emits chromosomes in reference (``.fai``) order, so concatenating the
   per-chromosome row files in ``.fai`` order reproduces the serial file. The same
   applies to the triplet counts: ``triplets.json`` preserves dict insertion
   order, so the order in which contexts are first seen is part of the output and
   the per-chromosome counts must be merged in ``.fai`` order too -- not in the
   order the pool happened to return them.

The guarantee is on the *decompressed* content of ``matching_bases.csv.gz``: the
exact bytes of a gzip member are not reproducible across gzip versions. The final
file is still written by a single writer in the parent process, exactly as before.
"""

from __future__ import annotations

import csv
import gzip
import json
import os
import shutil
import subprocess
import tempfile
from collections import defaultdict
from io import TextIOWrapper

from .mutation_extractor_manager import _chroms_with_reads, _read_fai_chroms
from .utils import log
from .parallel_mpileup import mpileup_cmd
from .parallel_resources import cap_jobs, get_mp_context

#: Rough resident memory per scan worker (mpileup child + the Python buffer).
#: Only used to cap the pool when the caller passes a memory budget.
SCAN_MB_PER_JOB = 1500


def scan_stream(extractor, stream, write_row):
    """Run the inherited 3-line sliding window over one pileup text stream.

    This is a faithful transcription of the loop in
    ``MultipleSpeciesMutationExtractor.extract``; every decision it makes is
    delegated to the (inherited) methods of ``extractor``. Returns the triplet
    context counts for this stream.
    """
    triplet_counts = defaultdict(int)

    buffer = [None,
              extractor._parse_line(stream.readline()),
              extractor._parse_line(stream.readline())]
    qc_flags = [False,
                extractor._quality_check(buffer[1]),
                extractor._quality_check(buffer[2])]

    for line in stream:
        buffer = [buffer[1], buffer[2], extractor._parse_line(line)]
        qc_flags = [qc_flags[1], qc_flags[2], extractor._quality_check(buffer[2])]
        if all(qc_flags) and extractor._consecutive(*buffer):
            result, triplet = extractor._detect_site(buffer)
            if triplet is not None:
                triplet_counts[triplet] += 1
            if result is not None:
                write_row(result)

    return triplet_counts


def _scan_region(chrom, ref_fasta, bams, n_species, species_list, mapping, rows_path):
    """Worker: pileup one chromosome from the indexed BAMs and scan it.

    Returns ``(chrom, triplet_counts, rows_path or None)``. The detected rows are
    written to ``rows_path`` as plain text rather than returned, so a worker's
    memory stays bounded no matter how many sites a chromosome carries; the
    parent concatenates those files into the single gzip stream.
    """
    # Imported here, not at module scope: the extractor module imports this one,
    # so a top-level import would be circular.
    from .multiple_species_mutation_extractor_manager import (
        MultipleSpeciesMutationExtractor,
    )

    # A scan-only extractor: tree=None means it does detection only, and
    # output_dir is the temp directory, so nothing of the real run is touched.
    extractor = MultipleSpeciesMutationExtractor(
        pileup_file=None,
        output_dir=os.path.dirname(rows_path),
        n_species=n_species,
        tree=None,
        species_list=species_list,
        mapping=mapping,
        no_cache=True,
        verbose=False,
    )

    # Same options as Pileup.generate, plus -r: see parallel_mpileup.
    cmd = mpileup_cmd(ref_fasta, bams, region=chrom)
    # stderr goes to a file, never to DEVNULL: a failing mpileup must be able to
    # say why, naming the chromosome. A file (rather than a second pipe) also
    # rules out the classic deadlock of draining one pipe while the other fills.
    stderr_path = rows_path + ".stderr"
    with open(stderr_path, "w+") as stderr_file:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=stderr_file)
        try:
            with open(rows_path, "w", newline="") as rows_file:
                writer = csv.writer(rows_file)
                with TextIOWrapper(proc.stdout) as stream:
                    triplet_counts = scan_stream(extractor, stream, writer.writerow)
        finally:
            proc.wait()
        stderr_file.seek(0)
        stderr_text = stderr_file.read()
    os.remove(stderr_path)

    if proc.returncode != 0:
        raise RuntimeError(
            f"samtools mpileup failed for chromosome {chrom} "
            f"(exit {proc.returncode}):\n{stderr_text.strip()}"
        )

    if os.path.getsize(rows_path) == 0:
        os.remove(rows_path)
        rows_path = None
    return chrom, dict(triplet_counts), rows_path


def scan_pileup_parallel(extractor, csv_path, triplets_path, header):
    """Scan the pileup one chromosome at a time and write both output files.

    ``extractor`` is the :class:`MultipleSpeciesMutationExtractor` driving the
    run; its ``ref_fasta`` / ``bams`` / ``fai_path`` / ``cores`` / ``scan_jobs``
    / ``max_memory_mb`` attributes configure the split, and its detection methods
    do the actual work inside each worker. Returns the merged triplet counts,
    having already written ``csv_path`` and ``triplets_path``.
    """
    chrom_lengths = _read_fai_chroms(extractor.fai_path)
    chrom_order = [c for c, _ in chrom_lengths]
    with_reads = _chroms_with_reads(extractor.bams)
    # Longest chromosome first so the big ones start early and dominate the
    # makespan; chromosomes without reads emit nothing and are skipped.
    tasks = [c for c, _ in sorted(chrom_lengths, key=lambda cl: cl[1], reverse=True)
             if c in with_reads]

    if not tasks:
        raise RuntimeError(
            f"No chromosome in {extractor.fai_path} carries mapped reads in any of "
            f"{len(extractor.bams)} BAM(s); nothing to scan."
        )

    n_workers = cap_jobs(extractor.scan_jobs, extractor.cores, len(tasks),
                         mb_per_job=SCAN_MB_PER_JOB,
                         max_memory_mb=extractor.max_memory_mb)
    log(f"Scanning pileup in parallel: {len(tasks)} chromosomes with reads "
        f"over {n_workers} workers...", extractor.verbose)

    tmp_dir = tempfile.mkdtemp(prefix="coral_scan_", dir=extractor.output_dir)
    try:
        args = [(c, extractor.ref_fasta, extractor.bams, extractor.n_species,
                 extractor.species_list, extractor.mapping,
                 os.path.join(tmp_dir, f"{i:05d}.rows"))
                for i, c in enumerate(tasks)]

        if n_workers <= 1:
            results = [_scan_region(*a) for a in args]
        else:
            # chunksize=1: hand out one chromosome per grab, so the few big arms
            # land on different workers. The default chunksize would bundle
            # contiguous (here longest-first) tasks onto one worker and serialise
            # the heavy work.
            with get_mp_context().Pool(n_workers) as pool:
                results = pool.starmap(_scan_region, args, chunksize=1)

        triplets_by_chrom = {}
        rows_by_chrom = {}
        for chrom, triplets, rows_path in results:
            triplets_by_chrom[chrom] = triplets
            if rows_path is not None:
                rows_by_chrom[chrom] = rows_path

        # Accumulate in .fai order, NOT in the order the pool returned results
        # (which is longest-chromosome-first). triplets.json preserves dict
        # insertion order, so the order contexts are FIRST SEEN is part of the
        # output -- merging in task order silently reorders the file.
        triplet_counts = defaultdict(int)
        for chrom in chrom_order:
            for context, count in triplets_by_chrom.get(chrom, {}).items():
                triplet_counts[context] += count

        # A single writer emitting .fai order reproduces the serial scan's
        # content exactly. Written to a temporary name and renamed, so a run
        # killed mid-merge cannot leave a truncated file that the cache check
        # would later accept as complete.
        tmp_csv = csv_path + ".tmp"
        with gzip.open(tmp_csv, "wt", newline="") as outfile:
            csv.writer(outfile).writerow(header)
            for chrom in chrom_order:
                rows_path = rows_by_chrom.get(chrom)
                if rows_path is None:
                    continue
                with open(rows_path, "r", newline="") as rows_file:
                    shutil.copyfileobj(rows_file, outfile)
        os.replace(tmp_csv, csv_path)
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)

    triplet_counts = dict(triplet_counts)
    with open(triplets_path, "w") as tf:
        json.dump(triplet_counts, tf, indent=2)
    log(f"Saved {len(triplet_counts)} triplet contexts to {triplets_path}",
        extractor.verbose)
    return triplet_counts
