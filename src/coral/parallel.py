"""Parallel execution utilities for the multi-species pipeline.

CPU budgets respect the process affinity mask when available, so parallelism
matches the resources granted by the host or scheduler.

All pools use the `spawn` start method for consistent behaviour across
platforms. Worker functions must therefore be module-level and use picklable
arguments, and any state needed by the parent must be returned explicitly.
"""

from __future__ import annotations

import csv
import gzip
import json
import multiprocessing
import os
import shutil
import subprocess
import tempfile
from collections import defaultdict
from io import TextIOWrapper

from .pileup_manager import mpileup_cmd
from .utils import log


# --------------------------------------------------------------------------
# Resources: CPU/memory budgets and the multiprocessing context
# --------------------------------------------------------------------------

#: Start method used for every pool in this module.
START_METHOD = "spawn"


def get_mp_context():
    """Return the multiprocessing context every pool in this module must use."""
    return multiprocessing.get_context(START_METHOD)


def detect_cores():
    """Number of CPUs this *process* may use, not the number the host owns."""
    if hasattr(os, "sched_getaffinity"):
        try:
            return max(1, len(os.sched_getaffinity(0)))
        except OSError:
            pass
    return max(1, multiprocessing.cpu_count())


def resolve_cores(cores=None):
    """Explicit ``--cores`` wins; otherwise auto-detect."""
    if cores is not None:
        if cores < 1:
            raise ValueError(f"cores must be >= 1, got {cores}")
        return int(cores)
    return detect_cores()


def cap_jobs(requested, cores, n_tasks, mb_per_job=None, max_memory_mb=None):
    """Workers for one stage: min(requested-or-cores, cores, tasks, memory ceiling).

    ``mb_per_job`` is that stage's own estimate of resident memory per worker.
    """
    n = cores if requested is None else int(requested)
    if n < 1:
        raise ValueError(f"job count must be >= 1, got {requested}")
    if n > cores:
        log(f"Requested {n} jobs but only {cores} core(s) available; using {cores}.")
        n = cores
    n = min(n, max(1, int(n_tasks)))
    if max_memory_mb and mb_per_job:
        memory_jobs = int(max_memory_mb // mb_per_job)
        if memory_jobs < 1:
            # mb_per_job is a rough estimate, so warn rather than fail the run.
            log(f"Memory budget {max_memory_mb} MB is below the ~{mb_per_job} MB "
                f"estimated per worker; running 1 worker anyway.")
        n = min(n, max(1, memory_jobs))
    return max(1, n)


# --------------------------------------------------------------------------
# Alignment: one worker per species
# --------------------------------------------------------------------------
# Each species is aligned independently. Although `bwa mem` is threaded, the
# downstream SAM filtering is single-threaded, so running species concurrently
# makes better use of the available cores.
#
# Approximate memory per alignment worker, used only when a memory limit is set.
ALIGN_MB_PER_JOB = 6000


def _init_align_worker():
    # The SAM filters draw MAPQ histograms; workers have no display.
    import matplotlib
    matplotlib.use("Agg")


def _align_one(aligner, streamed, kwargs):
    """Worker: run one species' alignment, return its BAM path and filter stats. """
    if streamed:
        aligner.align_streamed(**kwargs)
    else:
        aligner.align_disk_cached(**kwargs)
    return aligner.species, aligner.final_bam, aligner.filter_stats


def threads_per_job(cores, n_jobs):
    """Split the core budget across concurrent aligners, at least one each."""
    return max(1, cores // max(1, n_jobs))


def split_threads(cores, n_workers, n_tasks):
    """Per-task thread counts, with the division remainder handed to the first few tasks. """
    base = threads_per_job(cores, n_workers)
    if n_workers != n_tasks:
        return [base] * n_tasks
    extra = max(0, cores - base * n_workers)
    return [base + (1 if i < extra else 0) for i in range(n_tasks)]


def _input_size(aligner):
    """Fragment FASTQ size, as a stand-in for how long a species will take."""
    try:
        return os.path.getsize(aligner.fastq)
    except (OSError, TypeError):
        return 0


def plan_align_jobs(n_species, cores, align_jobs=None, max_memory_mb=None):
    """(workers, per-species thread counts) for the alignment stage."""
    n_workers = cap_jobs(align_jobs, cores, n_species,
                         mb_per_job=ALIGN_MB_PER_JOB, max_memory_mb=max_memory_mb)
    return n_workers, split_threads(cores, n_workers, n_species)


def run_alignments(aligners, streamed, align_kwargs, n_workers, verbose=True):
    """Align every species in a pool and return the aligners in *input* order. """
    if not aligners:
        return []

    payload = [(aligner, streamed, align_kwargs) for aligner in aligners]

    if n_workers <= 1:
        log("Aligning species sequentially...", verbose)
        results = [_align_one(*args) for args in payload]
    else:
        # Submit the biggest inputs first - the returned order is unchanged.
        payload.sort(key=lambda item: _input_size(item[0]), reverse=True)
        t = sorted({a.cores for a in aligners})
        spread = str(t[0]) if len(t) == 1 else f"{t[0]}-{t[-1]}"
        log(f"Aligning {len(aligners)} species over {n_workers} workers "
            f"({spread} thread(s) each)...", verbose)
        with get_mp_context().Pool(n_workers, initializer=_init_align_worker) as pool:
            results = pool.starmap(_align_one, payload, chunksize=1)

    produced = {species: (path, stats) for species, path, stats in results}
    for aligner in aligners:
        path, stats = produced.get(aligner.species, (None, None))
        if path is None or not os.path.exists(path):
            raise FileNotFoundError(
                f"Alignment for {aligner.species} did not produce {aligner.final_bam}"
            )
        # Carry the worker's diagnostics onto the parent's copy, so the serial
        # and parallel paths leave the aligners in the same state.
        aligner.filter_stats = stats

    return list(aligners)


# --------------------------------------------------------------------------
# Scan: one worker per chromosome
# --------------------------------------------------------------------------
# The serial scan reads a whole-genome pileup with a 3-line sliding window. Each
# worker here pileups one chromosome straight from the indexed BAMs and runs the
# same window over it, using the same mpileup options (see `mpileup_cmd`).


#: Rough resident memory per scan worker (mpileup child + the Python buffer).
#: Only used to cap the pool when the caller passes a memory budget.
SCAN_MB_PER_JOB = 1500


def scan_stream(extractor, stream, write_row):
    """Run the 3-line sliding window over one pileup text stream.

    A transcription of the loop in ``MultipleSpeciesMutationExtractor.extract``;
    every decision is delegated to that extractor's own methods. Returns this
    stream's triplet context counts.
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
    """Worker: pileup one chromosome and scan it.

    Returns ``(chrom, triplet_counts, rows_path or None)``. Detected rows go to
    ``rows_path`` rather than into the return value, so a worker's memory stays
    bounded however many sites a chromosome carries; the parent concatenates
    those files into the single gzip stream.
    """

    from .multiple_species_mutation_extractor_manager import (
        MultipleSpeciesMutationExtractor,
    )
    
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

    cmd = mpileup_cmd(ref_fasta, bams, region=chrom)
    # stderr to a file to say tracing the reason and chromosome failure.
    stderr_path = rows_path + ".stderr"
    with open(stderr_path, "w+") as stderr_file:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=stderr_file)
        try:
            with open(rows_path, "w", newline="") as rows_file:
                writer = csv.writer(rows_file)
                with TextIOWrapper(proc.stdout) as stream:
                    triplet_counts = scan_stream(extractor, stream, writer.writerow)
        except BaseException:
            # Don't wait on a child whose output we've stopped reading: a failed
            # job would otherwise leave samtools hanging on the cluster.
            if proc.stdout is not None:
                proc.stdout.close()
            proc.terminate()
            try:
                proc.wait(timeout=5)
            except subprocess.TimeoutExpired:
                proc.kill()
                proc.wait()
            raise
        else:
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

    ``extractor`` drives the run: its ``ref_fasta`` / ``bams`` / ``fai_path`` /
    ``cores`` / ``scan_jobs`` / ``max_memory_mb`` configure the split, and its
    detection methods do the work inside each worker. Returns the merged triplet
    counts, having written ``csv_path`` and ``triplets_path``.
    """
    # Imported here to keep this module light for the workers that don't need it.
    from .mutation_extractor_manager import _chroms_with_reads, _read_fai_chroms

    chrom_lengths = _read_fai_chroms(extractor.fai_path)
    chrom_order = [c for c, _ in chrom_lengths]
    with_reads = _chroms_with_reads(extractor.bams)
    # Longest first so the big ones start early and dominate the makespan;
    # chromosomes without reads emit nothing and are skipped.
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

    # One generation with the CSV: clear stale counts so a crash between the two
    # writes can't pair a fresh CSV with old triplets.
    if os.path.exists(triplets_path):
        os.remove(triplets_path)

    tmp_dir = tempfile.mkdtemp(prefix="coral_scan_", dir=extractor.output_dir)
    try:
        args = [(c, extractor.ref_fasta, extractor.bams, extractor.n_species,
                 extractor.species_list, extractor.mapping,
                 os.path.join(tmp_dir, f"{i:05d}.rows"))
                for i, c in enumerate(tasks)]

        if n_workers <= 1:
            results = [_scan_region(*a) for a in args]
        else:
            with get_mp_context().Pool(n_workers) as pool:
                results = pool.starmap(_scan_region, args, chunksize=1)

        triplets_by_chrom = {}
        rows_by_chrom = {}
        for chrom, triplets, rows_path in results:
            triplets_by_chrom[chrom] = triplets
            if rows_path is not None:
                rows_by_chrom[chrom] = rows_path

        # Accumulate in .fai order, NOT in the order the pool returned results
        # (longest-chromosome-first).
        triplet_counts = defaultdict(int)
        for chrom in chrom_order:
            for context, count in triplets_by_chrom.get(chrom, {}).items():
                triplet_counts[context] += count

        # A single writer emitting .fai order reproduces the serial scan's
        # content exactly. Written to a temporary name and renamed to avoid trucation.
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
    # Renamed like the CSV, so a kill mid-write can't leave truncated JSON that
    # the cache check would then accept.
    tmp_triplets = triplets_path + ".tmp"
    with open(tmp_triplets, "w") as tf:
        json.dump(triplet_counts, tf, indent=2)
    os.replace(tmp_triplets, triplets_path)
    log(f"Saved {len(triplet_counts)} triplet contexts to {triplets_path}",
        extractor.verbose)
    return triplet_counts


# --------------------------------------------------------------------------
# Fitch: flattened tree, optionally over row blocks
# --------------------------------------------------------------------------
# The serial path pays an ete3 deep copy of the tree *per site*, only because
# `_recursive_state_check` mutates the tree by attaching a `state` feature. Two
# output-preserving changes are made here:
#
# (a) Flatten the tree once into plain arrays. Fitch then writes states into a
#     per-row scratch list, so nothing is mutated and no copy is needed.
# (b) Parallelise over row blocks. 

#: Rows per block handed to a worker. 
FITCH_CHUNK_SIZE = 20_000

#: Rough resident memory per Fitch worker (one DataFrame block + its results).
FITCH_MB_PER_JOB = 800


class FlatTree:
    """A Fitch-ready, picklable, immutable view of an annotated CORAL tree.

    Node ids follow the same pre-order the original recursion visits.
    """

    __slots__ = ("children", "custom_name", "parent_name", "leaf_col", "root", "taxa_columns")

    def __init__(self, tree, mapping):
        self.children = []
        self.custom_name = []
        self.parent_name = []
        self.leaf_col = []      # index into the row's taxa values, or -1 for internal
        self.root = 0

        # The `taxaK` column order of matching_bases.csv.gz, as built by the
        # extractor: the integer keys of the mapping in insertion order.
        self.taxa_columns = [f"taxa{k}" for k in mapping if isinstance(k, int)]
        col_index = {name: i for i, name in enumerate(self.taxa_columns)}

        def add(node):
            idx = len(self.children)
            self.children.append([])
            self.custom_name.append(node.custom_name)
            self.parent_name.append(node.up.custom_name if node.up else "ROOT")
            if node.is_leaf():
                self.leaf_col.append(col_index[f"taxa{mapping[node.name]}"])
            else:
                self.leaf_col.append(-1)
            for child in node.children:
                self.children[idx].append(add(child))
            return idx

        add(tree)

    def __reduce__(self):
        # __slots__ without __dict__: define pickling explicitly for spawn.
        return (_rebuild_flat_tree, (self.children, self.custom_name,
                                     self.parent_name, self.leaf_col,
                                     self.root, self.taxa_columns))


def _rebuild_flat_tree(children, custom_name, parent_name, leaf_col, root, taxa_columns):
    flat = FlatTree.__new__(FlatTree)
    flat.children = children
    flat.custom_name = custom_name
    flat.parent_name = parent_name
    flat.leaf_col = leaf_col
    flat.root = root
    flat.taxa_columns = taxa_columns
    return flat


def _states(flat, taxa_values):
    """Fitch downward pass: the state set of every node for one site.

    Sets are never mutated in place, so a single-child node legitimately shares
    its child's set object, as the original does when it assigns
    ``node_state = child_states[0]``.
    """
    states = [None] * len(flat.children)

    def visit(i):
        col = flat.leaf_col[i]
        if col >= 0:
            states[i] = {taxa_values[col]}
            return states[i]
        child_states = [visit(c) for c in flat.children[i]]
        node_state = child_states[0]
        for child_state in child_states[1:]:
            intersect = node_state & child_state
            node_state = intersect if intersect else node_state | child_state
        states[i] = node_state
        return node_state

    visit(flat.root)
    return states


def _descend(flat, states, i, parent_state, chrom, pos, left, right, out):
    """Fitch upward pass for one site; returns this subtree's ambiguity count."""
    next_state = parent_state
    if parent_state not in states[i]:
        if len(states[i]) > 1:
            # Ambiguous: the original returns here, so the subtree is not walked.
            return 1
        next_state = next(iter(states[i]))
        branch_key = f"{flat.parent_name[i]}→{flat.custom_name[i]}"
        mutation = f"{left}[{parent_state}>{next_state}]{right}"
        out.setdefault(branch_key, []).append((chrom, pos, mutation))
    ambiguous = 0
    for child in flat.children[i]:
        ambiguous += _descend(flat, states, child, next_state, chrom, pos, left, right, out)
    return ambiguous


def fitch_block(flat, block):
    """Run Fitch over one row block. Returns ``(mutation_dict, ambiguous_count)``."""
    out = {}
    ambiguous = 0
    n_meta = 4  # chromosome, position, left, right
    columns = ["chromosome", "position", "left", "right"] + flat.taxa_columns
    for row in block[columns].itertuples(index=False, name=None):
        chrom, pos, left, right = row[:n_meta]
        states = _states(flat, row[n_meta:])
        root_state = states[flat.root]
        if len(root_state) != 1:
            ambiguous += 1
            continue
        ambiguous += _descend(flat, states, flat.root, next(iter(root_state)),
                              chrom, pos, left, right, out)
    return out, ambiguous


_WORKER_FLAT = None


def _init_fitch_worker(flat):
    global _WORKER_FLAT
    _WORKER_FLAT = flat


def _fitch_block_worker(block):
    return fitch_block(_WORKER_FLAT, block)


def run_fitch(csv_path, tree, mapping, cores=1, fitch_jobs=None, max_memory_mb=None,
              chunk_size=FITCH_CHUNK_SIZE, log_fn=None):
    """Flattened Fitch over ``matching_bases.csv.gz``, in one process or a pool.

    Returns ``(mutation_dict, ambiguous_count)`` ready for ``_save_results``.
    With ``fitch_jobs=1`` this is the flattened single-process implementation,
    still much faster than the per-row ``tree.copy()``; the untouched original
    recursion stays available through the extractor's serial path.
    """
    # Imported here so alignment and scan workers don't pay for pandas.
    import pandas as pd

    if chunk_size < 1:
        raise ValueError(f"chunk_size must be >= 1, got {chunk_size}")
    flat = FlatTree(tree, mapping)
    # Task count is unknown until the file is read; cap on cores and let the pool
    # size be the ceiling rather than the exact block count.
    n_workers = cap_jobs(fitch_jobs, cores, cores,
                         mb_per_job=FITCH_MB_PER_JOB, max_memory_mb=max_memory_mb)

    blocks = pd.read_csv(csv_path, chunksize=chunk_size)
    mutation_dict = {}
    ambiguous_total = 0

    def absorb(result):
        nonlocal ambiguous_total
        out, ambiguous = result
        ambiguous_total += ambiguous
        for branch_key, records in out.items():
            # imap preserves input order and blocks are contiguous, so appending
            # here reproduces the serial per-branch record order.
            mutation_dict.setdefault(branch_key, []).extend(records)

    if n_workers <= 1:
        if log_fn:
            log_fn("Running Fitch reconstruction in a single process (flattened tree)...")
        for block in blocks:
            absorb(fitch_block(flat, block))
    else:
        if log_fn:
            log_fn(f"Running Fitch reconstruction over {n_workers} workers...")
        # One pool for the whole file (not one per block) so a worker's spawn
        # cost is paid once; the flat tree ships via the initializer.
        with get_mp_context().Pool(n_workers, initializer=_init_fitch_worker, initargs=(flat,)) as pool:
            for result in pool.imap(_fitch_block_worker, blocks):
                absorb(result)

    return mutation_dict, ambiguous_total
