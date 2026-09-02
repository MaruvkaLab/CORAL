"""Portable resolution of CPU/memory budgets and of the multiprocessing context.

The same code has to behave correctly on a personal machine and inside an HTCondor
slot, so nothing here assumes that the whole host belongs to us:

* CPU count comes from the process affinity mask when the platform exposes one
  (Linux, and therefore Condor, where the slot's mask is what we were actually
  granted); ``multiprocessing.cpu_count()`` is only the fallback (macOS has no
  ``sched_getaffinity``).
* An optional memory budget caps the number of concurrent workers per stage,
  because parallelism multiplies resident memory and an over- or under-sized
  request is penalised on the cluster.
* The start method is pinned to ``spawn`` on every platform so that worker
  behaviour is identical on macOS (spawn by default) and Linux (fork by default),
  and so that no worker can silently depend on inherited global state.
"""

from __future__ import annotations

import multiprocessing
import os

#: Start method used for every pool in this package. Pinned rather than inherited
#: so macOS and Linux behave identically; all worker callables are therefore
#: module-level functions taking picklable arguments only.
START_METHOD = "spawn"


def get_mp_context():
    """Return the multiprocessing context every pool in this package must use."""
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
    """Explicit ``--cores`` wins; otherwise auto-detect.

    Unlike CORAL's single-species pipeline (where ``cores=None`` means serial),
    an unspecified core count here means "use what you were given" -- parallel is
    the default for ``run_multi`` and is switched off explicitly. See
    :class:`~coral.pipeline.MultiSpeciesMutationPipeline`.
    """
    if cores is not None:
        if cores < 1:
            raise ValueError(f"cores must be >= 1, got {cores}")
        return int(cores)
    return detect_cores()


def cap_jobs(requested, cores, n_tasks, mb_per_job=None, max_memory_mb=None):
    """Workers for one stage: min(requested-or-cores, tasks, memory ceiling).

    ``mb_per_job`` is that stage's own estimate of resident memory per worker; it
    is only consulted when the caller supplied ``max_memory_mb``. Returns at
    least 1, so a stage never ends up with an empty pool.
    """
    n = cores if requested is None else int(requested)
    if n < 1:
        raise ValueError(f"job count must be >= 1, got {requested}")
    n = min(n, max(1, int(n_tasks)))
    if max_memory_mb and mb_per_job:
        n = min(n, max(1, int(max_memory_mb // mb_per_job)))
    return max(1, n)
