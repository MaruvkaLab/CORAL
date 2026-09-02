"""Species-parallel driver for CORAL's existing :class:`~coral.alignment_manager.Aligner`.

The multi-species pipeline aligns species one after another. The aligner command
itself is threaded (``bwa mem -t``), but the SAM filter that consumes its output
-- ``with_continuity_filter_sam`` / ``filter_sam`` -- is a single-threaded Python
loop over every read, and it is what throttles the stage. Running ``S`` species
concurrently with ``cores // S`` threads each therefore wins: the alignment
threads are traded for S concurrent filter loops.

Nothing about a species' alignment depends on another species: ``final_bam``,
``raw_bam``, ``hist_name`` and ``log_path`` are all per-species (see
``Aligner.__init__``), and the shared ``BAMs``/``Plots`` directories are only ever
created with ``exist_ok=True``. The per-species outputs are therefore identical
to the serial run.

The one ordering hazard: a pool returns results as they finish, but the BAM order
is what defines the ``taxaK`` columns downstream. This module always returns the
aligners in the caller's input order, never in completion order.
"""

from __future__ import annotations

import os

from .utils import log

from .parallel_resources import cap_jobs, get_mp_context

#: Rough resident memory per alignment worker: the aligner holds the reference
#: index in RAM, which dominates. Only consulted when a memory budget is given.
ALIGN_MB_PER_JOB = 6000


def _init_worker():
    # The SAM filters draw MAPQ histograms; workers have no display.
    import matplotlib
    matplotlib.use("Agg")


def _align_one(aligner, streamed, kwargs):
    """Worker: run one species' alignment and return its BAM path.

    ``aligner`` is a plain data object (paths, ints, a command template), so it
    round-trips through pickle unchanged; the parent keeps its own copy, whose
    ``final_bam`` is the same path this worker writes.
    """
    if streamed:
        aligner.align_streamed(**kwargs)
    else:
        aligner.align_disk_cached(**kwargs)
    return aligner.species, aligner.final_bam


def threads_per_job(cores, n_jobs):
    """Split the core budget across concurrent aligners, at least one each."""
    return max(1, cores // max(1, n_jobs))


def plan_jobs(n_species, cores, align_jobs=None, max_memory_mb=None):
    """(workers, threads per worker) for the alignment stage."""
    n_workers = cap_jobs(align_jobs, cores, n_species,
                         mb_per_job=ALIGN_MB_PER_JOB, max_memory_mb=max_memory_mb)
    return n_workers, threads_per_job(cores, n_workers)


def run_alignments(aligners, streamed, align_kwargs, n_workers, verbose=True):
    """Align every species, in a pool, and return the aligners in input order.

    ``aligners`` must already be built by the caller (same construction as the
    serial pipeline, but with the reduced per-worker thread count).
    """
    if not aligners:
        return []

    payload = [(aligner, streamed, align_kwargs) for aligner in aligners]

    if n_workers <= 1:
        log("Aligning species sequentially...", verbose)
        results = [_align_one(*args) for args in payload]
    else:
        log(f"Aligning {len(aligners)} species over {n_workers} workers "
            f"({aligners[0].cores} thread(s) each)...", verbose)
        with get_mp_context().Pool(n_workers, initializer=_init_worker) as pool:
            # starmap returns results in input order, but we do not rely on that:
            # the aligner list below is the caller's, so column order is anchored
            # to species_list regardless of completion order.
            results = pool.starmap(_align_one, payload, chunksize=1)

    produced = dict(results)
    for aligner in aligners:
        path = produced.get(aligner.species)
        if path is None or not os.path.exists(path):
            raise FileNotFoundError(
                f"Alignment for {aligner.species} did not produce {aligner.final_bam}"
            )

    return list(aligners)
