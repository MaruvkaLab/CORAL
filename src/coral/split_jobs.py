"""``coral run_multi --split-jobs``: run the multi-species pipeline as an HTCondor
DAG, with one job per ingroup species for the alignment. Without --split-jobs
nothing in this module runs.

Every DAG node runs the same ``coral run_multi`` command plus a hidden
``--split-stage``:

    prepare      download and index the outgroup
    align_NNNN   download, fragment and align ONE ingroup species
    main         the ordinary run_multi; every alignment is cached, so it goes
                 straight on to the scan, PHYLIP and Fitch

The align stage makes the same Genome/Aligner calls as the in-process pipeline;
only the thread count differs, which does not change single-end alignments.
No file-transfer commands are used: the output directory must be shared.
"""

import os
import shlex
import shutil
import subprocess
import sys

import pysam

from .alignment_manager import Aligner
from .genome_manager import Genome
from .run_phylip import check_phylip_available
from .utils import log

#: Copied into the jobs, which do not inherit the submitter's environment.
#: Secrets such as NCBI_API_KEY are not: run_stage.sh is written into the run dir.
PASSED_ENV = [
    "PATH", "LD_LIBRARY_PATH", "CONDA_PREFIX", "HOME", "MPLBACKEND",
    "OMP_NUM_THREADS", "OMP_THREAD_LIMIT", "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "JAVA_TOOL_OPTIONS",
]

#: Index files Genome.index writes, per aligner (minimap2/bbmap index in memory).
INDEX_EXTS = {
    "bwa": ["amb", "ann", "bwt", "pac", "sa"],
    "bwa-mem2": ["0123", "amb", "ann", "bwt.2bit.64", "pac"],
}

#: Aligners that write one primary record per read (see _check_raw_bam).
ONE_RECORD_PER_READ = {"bwa", "bwa-mem2", "minimap2"}


def dispatch(pipeline, args, argv):
    """Entry point from the CLI when --split-jobs or --split-stage is given."""
    if args.split_stage == "main":
        pipeline.run()
    elif args.split_stage == "prepare":
        run_prepare_stage(pipeline)
    elif args.split_stage == "align":
        if args.split_index is None:
            raise ValueError("--split-stage align needs --split-index")
        run_align_stage(pipeline, args.split_index, align_threads(args.job_cores))
    else:
        submit(pipeline, args, argv)


# --- stages (run inside the cluster jobs) ----------------------------------

def align_threads(job_cores):
    """One CPU of the job is left for the samtools/filter processes beside the aligner."""
    if job_cores < 1:
        raise ValueError(f"--job-cores must be >= 1, got {job_cores}")
    return max(1, job_cores - 1)


def parse_species(pipeline):
    """Fill species_list / species_dict / outgroup_name exactly as run() does."""
    if pipeline.newick_tree:
        pipeline.parse_and_annotate_tree()
    else:
        pipeline.parse_and_annotate_list()
    return [s for s in pipeline.species_list if s != pipeline.outgroup_name]


def _genome(pipeline, species):
    # Same construction as download_index_and_fragment.
    return Genome(name=species, accession=pipeline.species_dict[species],
                  output_dir=pipeline.output_dir, no_cache=pipeline.no_cache,
                  verbose=pipeline.verbose)


def run_prepare_stage(pipeline):
    parse_species(pipeline)
    genome = _genome(pipeline, pipeline.outgroup_name)
    genome.download()
    genome.index(aligner=pipeline.aligner_name)
    log(f"Split stage 'prepare' done: {pipeline.outgroup_name} is indexed.", pipeline.verbose)


def _fragment(pipeline, genome):
    """generate_fragment_fastq as download_index_and_fragment calls it, but written
    under a temporary name, so a killed job leaves no partial FASTQ to be cached."""
    length = pipeline.params.get("fragment_length", 150)
    offset = pipeline.params.get("fragment_offset", 75)
    final = os.path.join(genome.output_dir, f"{genome.name}.fastq")
    if os.path.exists(final):
        return genome.generate_fragment_fastq(length=length, offset=offset)
    genome.generate_fragment_fastq(length=length, output_fastq=final + ".tmp",
                                   offset=offset, force=True)
    os.replace(final + ".tmp", final)
    genome.fastq_path = final


def run_align_stage(pipeline, index, threads):
    """Download, fragment and align ingroup species number ``index``."""
    species = parse_species(pipeline)
    if not 0 <= index < len(species):
        raise IndexError(f"--split-index {index} is out of range: the run has "
                         f"{len(species)} ingroup species")
    log(f"Split stage 'align' {index}: {species[index]}", pipeline.verbose)

    # The outgroup must be indexed already; otherwise every align job would
    # rebuild the same index files at once.
    reference = _genome(pipeline, pipeline.outgroup_name)
    needed = [reference.fasta_path + ext for ext in
              ["", ".fai"] + ["." + e for e in INDEX_EXTS.get((pipeline.aligner_name or "").lower(), [])]]
    missing = [p for p in needed if not os.path.exists(p)]
    if missing:
        raise FileNotFoundError("The outgroup is not prepared (run the 'prepare' stage "
                                "first); missing: " + ", ".join(missing))

    genome = _genome(pipeline, species[index])
    genome.download()
    _fragment(pipeline, genome)

    # Same construction and arguments as the alignment stage, except cores.
    aligner = Aligner(species_genome=genome, reference_genome=reference,
                      base_output_dir=pipeline.output_dir, aligner_cmd=pipeline.aligner_cmd,
                      aligner_name=pipeline.aligner_name, cores=threads,
                      verbose=pipeline.verbose)
    raw_was_cached = os.path.exists(aligner.raw_bam)
    kwargs = {"mapq": pipeline.params.get("mapq", 60),
              "low_mapq": pipeline.params.get("low_mapq", 1),
              "continuity": pipeline.params.get("continuity", True)}
    if pipeline.params.get("streamed", False):
        aligner.align_streamed(max_sort_mem=pipeline.params.get("max_samtools_mem", None), **kwargs)
    else:
        aligner.align_disk_cached(**kwargs)
        if not raw_was_cached and pipeline.aligner_cmd is None \
                and (pipeline.aligner_name or "").lower() in ONE_RECORD_PER_READ:
            _check_raw_bam(aligner, threads)
    _check_final_bam(aligner)
    log(f"Split stage 'align' {index} done: {aligner.final_bam}", pipeline.verbose)
    return aligner


def _set_aside(aligner, reason, raw=False):
    """Rename broken outputs to *.incomplete (never delete them) and fail."""
    paths = ([aligner.raw_bam] if raw else []) + [aligner.final_bam + s for s in ("", ".bai", ".csi")]
    moved = []
    for path in paths:
        if os.path.exists(path):
            os.replace(path, path + ".incomplete")
            moved.append(path + ".incomplete")
    raise RuntimeError(f"{reason}. Set aside so a rerun recomputes them: " + ", ".join(moved))


def _check_raw_bam(aligner, threads):
    """``aligner | samtools view`` runs without pipefail, so a killed aligner
    leaves a valid but short BAM. Every read must have one primary record."""
    records = int(pysam.view("-c", "-F", "0x900", "-@", str(threads), aligner.raw_bam))
    with open(aligner.fastq, "rb") as f:
        reads = sum(block.count(b"\n") for block in iter(lambda: f.read(1 << 24), b"")) // 4
    if records != reads:
        _set_aside(aligner, f"{aligner.raw_bam} has {records} primary records for {reads} "
                            f"reads: the aligner did not finish", raw=True)


def _check_final_bam(aligner):
    """A failed ``samtools sort`` still has its output renamed into place."""
    try:
        pysam.quickcheck(aligner.final_bam)
    except pysam.SamtoolsError:
        _set_aside(aligner, f"{aligner.final_bam} fails samtools quickcheck")


# --- submission (where condor_submit_dag works, e.g. a login node) ---------

def stage_script(argv, cwd, python, env, main_extra_args):
    """run_stage.sh, which every node runs: prepare | align INDEX | main."""
    run = f'exec {shlex.quote(python)} -m coral "${{ARGS[@]}}" --split-stage'
    return "\n".join(
        ["#!/bin/bash",
         "# Written by coral run_multi --split-jobs: run_stage.sh prepare | align INDEX | main",
         "set -euo pipefail"]
        + [f"export {name}={shlex.quote(value)}" for name, value in env]
        + [f"cd {shlex.quote(cwd)}",
           f"ARGS=({shlex.join(argv)})",
           'case "${1:-}" in',
           f"    prepare) {run} prepare ;;",
           f'    align)   {run} align --split-index "$2" ;;',
           f"    main)    {run} main {shlex.join(main_extra_args)} ;;",
           '    *) echo "usage: $0 prepare | align INDEX | main" >&2; exit 2 ;;',
           "esac", ""])


def submit_description(script, stage_args, cwd, log_path, cpus, memory_mb, extra):
    """No file-transfer commands. ``extra`` lines come last, so they win."""
    quoted = " ".join("'" + a.replace("'", "''") + "'" for a in [script] + stage_args)
    return "\n".join(
        ["universe       = vanilla",
         "executable     = /bin/bash",
         f'arguments      = "{quoted}"',
         f"initialdir     = {cwd}"]
        + [f"{kind:<14} = {log_path}.{ext}" for kind, ext in
           [("log", "log"), ("output", "out"), ("error", "err")]]
        + [f"request_cpus   = {cpus}", f"request_memory = {memory_mb}"]
        + list(extra or []) + ["queue", ""])


def dag_text(subs, species, max_concurrent):
    nodes = [f"align_{i:04d}" for i in range(len(species))]
    lines = [f"JOB prepare {subs['prepare']}"]
    for i, (node, name) in enumerate(zip(nodes, species)):
        lines += [f"# {node} = {name}", f"JOB {node} {subs['align']}",
                  f'VARS {node} index="{i}"', f"CATEGORY {node} align"]
    lines += [f"JOB main {subs['main']}",
              "PARENT prepare CHILD " + " ".join(nodes),
              "PARENT " + " ".join(nodes) + " CHILD main"]
    if max_concurrent:
        lines.append(f"MAXJOBS align {max_concurrent}")
    return "\n".join(lines + [""])


def submit(pipeline, args, argv):
    """Write the DAG for this run and submit it (unless --no-submit)."""
    if pipeline.no_cache:
        raise ValueError("--split-jobs cannot be combined with --no-cache: the "
                         "final job relies on the cached alignments.")
    if not check_phylip_available("dnapars"):
        raise RuntimeError("PHYLIP's dnapars, which the final job needs, is not on PATH.")
    if not args.no_submit and shutil.which("condor_submit_dag") is None:
        raise RuntimeError("condor_submit_dag not found (use --no-submit to only write the files).")
    threads = align_threads(args.job_cores)
    species = parse_species(pipeline)
    if not species:
        raise ValueError("No ingroup species to align.")

    job_dir = os.path.join(os.path.abspath(pipeline.output_dir), "split_jobs")
    log_dir = os.path.abspath(args.job_log_dir or os.path.join(job_dir, "logs"))
    os.makedirs(log_dir, exist_ok=True)
    os.makedirs(job_dir, exist_ok=True)
    dag = os.path.join(job_dir, f"{pipeline.run_id}.dag")
    if os.path.exists(dag + ".lock"):
        raise RuntimeError(f"{dag}.lock exists: a DAGMan for this run may still be running. "
                           "If condor_q shows none, the lock is stale; remove it and rerun.")

    env = [(name, os.environ[name]) for name in PASSED_ENV if os.environ.get(name)]
    pythonpath = [os.path.dirname(os.path.dirname(os.path.abspath(__file__)))]
    pythonpath += [p for p in os.environ.get("PYTHONPATH", "").split(os.pathsep) if p]
    env.append(("PYTHONPATH", os.pathsep.join(dict.fromkeys(pythonpath))))

    # The scan runs an mpileup beside each worker, hence 2 x --cores. Without
    # --cores the main job would take every core of its node, so 4 is passed.
    main_args = [] if args.cores else ["--cores", "4"]
    main_cpus = args.main_job_cores or 2 * (args.cores or 4)
    main_memory = args.main_job_memory_mb or args.max_memory_mb or 16384

    script = os.path.join(job_dir, "run_stage.sh")
    with open(script, "w") as f:
        f.write(stage_script(argv, os.getcwd(), sys.executable, env, main_args))
    os.chmod(script, 0o755)

    subs = {}
    for name, stage_args, log_name, cpus, memory, extra in [
            ("prepare", ["prepare"], "prepare", 1, args.job_memory_mb, args.job_extra),
            ("align", ["align", "$(index)"], "align_$(index)", args.job_cores,
             args.job_memory_mb, args.job_extra),
            ("main", ["main"], "main", main_cpus, main_memory, args.main_job_extra)]:
        subs[name] = os.path.join(job_dir, f"{name}.sub")
        with open(subs[name], "w") as f:
            f.write(submit_description(script, stage_args, os.getcwd(),
                                       os.path.join(log_dir, log_name), cpus, memory, extra))
    with open(dag, "w") as f:
        f.write(dag_text(subs, species, args.job_max_concurrent))

    log(f"Split run: 1 prepare job, {len(species)} align jobs ({args.job_cores} CPUs, "
        f"{threads} aligner threads, {args.job_memory_mb} MB each), 1 main job "
        f"({main_cpus} CPUs, {main_memory} MB). DAG: {dag}. Job logs: {log_dir}", pipeline.verbose)

    # -f rewrites DAGMan's files from an earlier submission; finished work is
    # not redone, because every stage reuses the cached files.
    cmd = ["condor_submit_dag", "-f", dag]
    if args.no_submit:
        log("--no-submit: nothing submitted. To submit:\n  " + shlex.join(cmd), pipeline.verbose)
    else:
        subprocess.run(cmd, cwd=job_dir, check=True)
    return dag
