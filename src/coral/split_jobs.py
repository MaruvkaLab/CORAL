"""``coral run_multi --split-jobs``: run the multi-species pipeline as an
HTCondor DAG, with one cluster job per ingroup species for the alignment.

Without ``--split-jobs`` nothing in this module runs.

The DAG has three kinds of node, each one the same ``coral run_multi`` command
with a hidden ``--split-stage`` appended:

    prepare        parse the tree/list, download and index the outgroup
    align_NNNN     download, fragment and align ONE ingroup species
    main           the ordinary ``run_multi``: every genome, FASTQ and BAM is
                   already cached, so it goes straight on to the scan, PHYLIP
                   and Fitch, through exactly the code a run without
                   ``--split-jobs`` uses

The align stage calls the same ``Genome`` and ``Aligner`` methods, with the same
arguments, as ``download_index_and_fragment`` and the alignment stage do, so
the files it leaves are the ones ``main`` would otherwise have written itself.
Only the thread count differs, and it does not change the alignments of
single-end fragment reads.

Only HTCondor is supported. The submit files use no file-transfer commands:
the output directory must be on a filesystem the execute nodes share.
"""

import os
import shlex
import shutil
import subprocess
import sys

from .alignment_manager import Aligner
from .genome_manager import Genome
from .run_phylip import check_phylip_available
from .utils import log

#: Environment variables copied into run_stage.sh when they are set at submit
#: time. Vanilla-universe jobs do not inherit the submitter's environment.
#: Secrets such as NCBI_API_KEY are deliberately not copied: the script is
#: written into the run directory.
PASSED_ENV = [
    "PATH", "LD_LIBRARY_PATH", "CONDA_PREFIX", "HOME", "MPLBACKEND",
    "OMP_NUM_THREADS", "OMP_THREAD_LIMIT", "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "JAVA_TOOL_OPTIONS",
]

#: Index files the align stage needs from the prepare stage, per aligner. These
#: are the files ``Genome.index`` checks; minimap2 and bbmap index in memory.
INDEX_EXTS = {
    "bwa": ["amb", "ann", "bwt", "pac", "sa"],
    "bwa-mem2": ["0123", "amb", "ann", "bwt.2bit.64", "pac"],
}

#: Aligners that write exactly one primary record per input read, so the raw
#: BAM can be checked against the FASTQ (see ``_check_raw_bam``).
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
    elif args.split_stage is None:
        submit(pipeline, args, argv)
    else:
        raise ValueError(f"Unknown --split-stage {args.split_stage!r}")


# --------------------------------------------------------------------------
# Stages (run inside the cluster jobs)
# --------------------------------------------------------------------------

def align_threads(job_cores):
    """Aligner/sort threads for a job with ``job_cores`` CPUs: one is left for
    the ``samtools view`` / SAM filter that run beside them."""
    if job_cores < 1:
        raise ValueError(f"--job-cores must be >= 1, got {job_cores}")
    return max(1, job_cores - 1)


def parse_species(pipeline):
    """Fill species_list / species_dict / outgroup_name exactly as run() does."""
    if pipeline.newick_tree:
        pipeline.parse_and_annotate_tree()
    else:
        pipeline.parse_and_annotate_list()


def ingroup(pipeline):
    """Ingroup species in species_list order; align job N aligns entry N."""
    return [s for s in pipeline.species_list if s != pipeline.outgroup_name]


def _genome(pipeline, species):
    # Same construction as download_index_and_fragment.
    return Genome(
        name=species,
        accession=pipeline.species_dict[species],
        output_dir=pipeline.output_dir,
        no_cache=pipeline.no_cache,
        verbose=pipeline.verbose,
    )


def run_prepare_stage(pipeline):
    """Download and index the outgroup, which every align job reads."""
    parse_species(pipeline)
    genome = _genome(pipeline, pipeline.outgroup_name)
    genome.download()
    genome.index(aligner=pipeline.aligner_name)
    log(f"Split stage 'prepare' done: outgroup {pipeline.outgroup_name} is "
        f"downloaded and indexed.", pipeline.verbose)


def _require_outgroup_index(reference, aligner_name):
    """The outgroup must already be indexed by the prepare job. Indexing it here
    would have every align job rebuild the same files at once."""
    exts = INDEX_EXTS.get((aligner_name or "").lower(), [])
    needed = [reference.fasta_path, f"{reference.fasta_path}.fai"]
    needed += [f"{reference.fasta_path}.{ext}" for ext in exts]
    missing = [p for p in needed if not os.path.exists(p)]
    if missing:
        raise FileNotFoundError(
            "The outgroup is not prepared (run the 'prepare' stage first); missing: "
            + ", ".join(missing))


def _fragment(pipeline, genome):
    """generate_fragment_fastq with the arguments download_index_and_fragment
    uses, but written to a temporary name and renamed when complete, so a job
    killed half-way does not leave a partial FASTQ that later runs would take
    as cached."""
    length = pipeline.params.get("fragment_length", 150)
    offset = pipeline.params.get("fragment_offset", 75)
    final = os.path.join(genome.output_dir, f"{genome.name}.fastq")
    if os.path.exists(final) and not pipeline.no_cache:
        return genome.generate_fragment_fastq(length=length, offset=offset,
                                              force=pipeline.no_cache)
    tmp = final + ".tmp"
    genome.generate_fragment_fastq(length=length, output_fastq=tmp, offset=offset,
                                   force=True)
    os.replace(tmp, final)
    genome.fastq_path = final
    return final


def run_align_stage(pipeline, index, threads):
    """Download, fragment and align ingroup species number ``index``."""
    parse_species(pipeline)
    species_names = ingroup(pipeline)
    if not 0 <= index < len(species_names):
        raise IndexError(f"--split-index {index} is out of range: the run has "
                         f"{len(species_names)} ingroup species")
    species = species_names[index]
    log(f"Split stage 'align' {index}: {species}", pipeline.verbose)

    reference = _genome(pipeline, pipeline.outgroup_name)
    _require_outgroup_index(reference, pipeline.aligner_name)

    genome = _genome(pipeline, species)
    genome.download()
    _fragment(pipeline, genome)

    # Same construction and arguments as the alignment stage, except cores.
    aligner = Aligner(
        species_genome=genome,
        reference_genome=reference,
        base_output_dir=pipeline.output_dir,
        aligner_cmd=pipeline.aligner_cmd,
        aligner_name=pipeline.aligner_name,
        cores=threads,
        verbose=pipeline.verbose,
    )
    raw_was_cached = os.path.exists(aligner.raw_bam)
    kwargs = {
        "mapq": pipeline.params.get("mapq", 60),
        "low_mapq": pipeline.params.get("low_mapq", 1),
        "continuity": pipeline.params.get("continuity", True),
    }
    if pipeline.params.get("streamed", False):
        kwargs["max_sort_mem"] = pipeline.params.get("max_samtools_mem", None)
        aligner.align_streamed(**kwargs)
    else:
        aligner.align_disk_cached(**kwargs)
        if not raw_was_cached and pipeline.aligner_cmd is None \
                and (pipeline.aligner_name or "").lower() in ONE_RECORD_PER_READ:
            _check_raw_bam(aligner, threads)

    _check_final_bam(aligner)
    log(f"Split stage 'align' {index} done: {aligner.final_bam}", pipeline.verbose)
    return aligner


def _set_aside(paths, reason):
    """Rename broken outputs to *.incomplete (never delete them) and fail."""
    moved = []
    for path in paths:
        if os.path.exists(path):
            os.replace(path, path + ".incomplete")
            moved.append(path + ".incomplete")
    raise RuntimeError(f"{reason}. Set aside so a rerun recomputes them: "
                       + ", ".join(moved))


def _fastq_records(path):
    lines = 0
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 24), b""):
            lines += block.count(b"\n")
    return lines // 4


def _check_raw_bam(aligner, threads):
    """The raw BAM is written by ``aligner | samtools view`` under /bin/sh, so
    an aligner killed mid-run (e.g. out of memory) still leaves a valid but
    short BAM that would be cached. Every read must have one primary record."""
    out = subprocess.run(
        ["samtools", "view", "-c", "-F", "0x900", "-@", str(threads), aligner.raw_bam],
        check=True, capture_output=True, text=True).stdout
    records = int(out.strip())
    reads = _fastq_records(aligner.fastq)
    if records != reads:
        _set_aside([aligner.raw_bam, aligner.final_bam,
                    aligner.final_bam + ".bai", aligner.final_bam + ".csi"],
                   f"{aligner.raw_bam} has {records} primary records for {reads} "
                   f"reads in {aligner.fastq}: the aligner did not finish")


def _check_final_bam(aligner):
    """A failed ``samtools sort`` still has its output renamed into place."""
    result = subprocess.run(["samtools", "quickcheck", aligner.final_bam])
    if result.returncode != 0:
        _set_aside([aligner.final_bam, aligner.final_bam + ".bai",
                    aligner.final_bam + ".csi"],
                   f"{aligner.final_bam} fails samtools quickcheck")


# --------------------------------------------------------------------------
# Submission (runs where condor_submit_dag works, e.g. a login node)
# --------------------------------------------------------------------------

def package_root():
    """Directory to put on PYTHONPATH so the jobs import this same coral."""
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def condor_quote(args):
    """HTCondor 'new syntax' arguments string."""
    return '"' + " ".join("'" + a.replace("'", "''") + "'" for a in args) + '"'


def stage_script(argv, cwd, python, env, main_extra_args):
    """Bash wrapper every node runs: run_stage.sh prepare | align N | main."""
    lines = [
        "#!/bin/bash",
        "# Written by coral run_multi --split-jobs. Runs one stage of the split run:",
        "#   run_stage.sh prepare | run_stage.sh align INDEX | run_stage.sh main",
        "set -euo pipefail",
    ]
    for name, value in env:
        lines.append(f"export {name}={shlex.quote(value)}")
    lines += [
        f"cd {shlex.quote(cwd)}",
        "ARGS=(" + " ".join(shlex.quote(a) for a in argv) + ")",
        'case "${1:-}" in',
        f'    prepare) exec {shlex.quote(python)} -m coral "${{ARGS[@]}}" --split-stage prepare ;;',
        f'    align)   exec {shlex.quote(python)} -m coral "${{ARGS[@]}}" --split-stage align --split-index "$2" ;;',
        f'    main)    exec {shlex.quote(python)} -m coral "${{ARGS[@]}}" --split-stage main'
        + "".join(" " + shlex.quote(a) for a in main_extra_args) + " ;;",
        '    *) echo "usage: $0 prepare | align INDEX | main" >&2; exit 2 ;;',
        "esac",
        "",
    ]
    return "\n".join(lines)


def submit_description(script, stage_args, cwd, log_dir, log_name, cpus, memory_mb,
                       extra):
    """One submit file. No file-transfer commands: the run dir must be shared.
    ``extra`` lines come last, so they override the defaults above them."""
    lines = [
        "# Written by coral run_multi --split-jobs.",
        "universe       = vanilla",
        "executable     = /bin/bash",
        f"arguments      = {condor_quote([script] + stage_args)}",
        f"initialdir     = {cwd}",
        f"log            = {os.path.join(log_dir, log_name)}.log",
        f"output         = {os.path.join(log_dir, log_name)}.out",
        f"error          = {os.path.join(log_dir, log_name)}.err",
        f"request_cpus   = {cpus}",
        f"request_memory = {memory_mb}",
    ]
    lines += list(extra or [])
    lines += ["queue", ""]
    return "\n".join(lines)


def dag_text(prepare_sub, align_sub, main_sub, species_names, max_concurrent):
    lines = ["# Written by coral run_multi --split-jobs.",
             f"JOB prepare {prepare_sub}"]
    names = []
    for i, species in enumerate(species_names):
        node = f"align_{i:04d}"
        names.append(node)
        lines += [f"# {node} = {species}",
                  f"JOB {node} {align_sub}",
                  f'VARS {node} index="{i}"',
                  f"CATEGORY {node} align"]
    lines.append(f"JOB main {main_sub}")
    lines.append("PARENT prepare CHILD " + " ".join(names))
    lines.append("PARENT " + " ".join(names) + " CHILD main")
    if max_concurrent:
        lines.append(f"MAXJOBS align {max_concurrent}")
    lines.append("")
    return "\n".join(lines)


def main_job_resources(args):
    """(request_cpus, request_memory, extra args) of the main job. The scan runs
    an mpileup beside each Python worker, so it asks for twice --cores. Without
    --cores the job would take every core on the node, so 4 is passed."""
    cores = args.cores
    extra_args = []
    if cores is None:
        cores = 4
        extra_args = ["--cores", str(cores)]
    cpus = args.main_job_cores or 2 * cores
    memory = args.main_job_memory_mb or args.max_memory_mb or 16384
    return cpus, memory, extra_args


def submit(pipeline, args, argv):
    """Write the DAG for this run and submit it (unless --no-submit)."""
    if pipeline.no_cache:
        raise ValueError("--split-jobs cannot be combined with --no-cache: the "
                         "final job relies on the cached alignments.")
    if not check_phylip_available("dnapars"):
        raise RuntimeError("PHYLIP's dnapars is required by the final job but was "
                           "not found on PATH.")
    if not args.no_submit and shutil.which("condor_submit_dag") is None:
        raise RuntimeError("condor_submit_dag was not found: --split-jobs needs "
                           "HTCondor (use --no-submit to only write the files).")
    threads = align_threads(args.job_cores)

    parse_species(pipeline)
    species_names = ingroup(pipeline)
    if not species_names:
        raise ValueError("No ingroup species to align.")

    out_dir = os.path.abspath(pipeline.output_dir)
    job_dir = os.path.join(out_dir, "split_jobs")
    log_dir = os.path.abspath(args.job_log_dir or os.path.join(job_dir, "logs"))
    os.makedirs(job_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)

    dag = os.path.join(job_dir, f"{pipeline.run_id}.dag")
    if os.path.exists(dag + ".lock"):
        raise RuntimeError(
            f"{dag}.lock exists: a DAGMan for this run may still be running. If "
            "condor_q shows none, the lock is stale; remove it and rerun.")

    cwd = os.getcwd()
    env = []
    for name in PASSED_ENV:
        if os.environ.get(name):
            env.append((name, os.environ[name]))
    pythonpath = [package_root()] + [p for p in os.environ.get("PYTHONPATH", "").split(os.pathsep) if p]
    env.append(("PYTHONPATH", os.pathsep.join(dict.fromkeys(pythonpath))))

    main_cpus, main_memory, main_extra_args = main_job_resources(args)

    script = os.path.join(job_dir, "run_stage.sh")
    with open(script, "w") as f:
        f.write(stage_script(argv, cwd, sys.executable, env, main_extra_args))
    os.chmod(script, 0o755)

    subs = {}
    for name, stage_args, log_name, cpus, memory, extra in [
        ("prepare", ["prepare"], "prepare", 1, args.job_memory_mb, args.job_extra),
        ("align", ["align", "$(index)"], "align_$(index)", args.job_cores,
         args.job_memory_mb, args.job_extra),
        ("main", ["main"], "main", main_cpus, main_memory, args.main_job_extra),
    ]:
        subs[name] = os.path.join(job_dir, f"{name}.sub")
        with open(subs[name], "w") as f:
            f.write(submit_description(script, stage_args, cwd, log_dir, log_name,
                                       cpus, memory, extra))

    with open(dag, "w") as f:
        f.write(dag_text(subs["prepare"], subs["align"], subs["main"], species_names,
                         args.job_max_concurrent))

    log(f"Split run: 1 prepare job, {len(species_names)} align jobs "
        f"({args.job_cores} CPUs, {threads} aligner threads, {args.job_memory_mb} MB each), "
        f"1 main job ({main_cpus} CPUs, {main_memory} MB).", pipeline.verbose)
    log(f"DAG: {dag}", pipeline.verbose)
    log(f"Job logs: {log_dir}", pipeline.verbose)

    # -f: rewrite DAGMan's own files from an earlier submission of this run.
    # Finished work is not redone: every stage reuses the cached files.
    cmd = ["condor_submit_dag", "-f", dag]
    if args.no_submit:
        log("--no-submit: nothing submitted. To submit:\n  " + " ".join(shlex.quote(c) for c in cmd),
            pipeline.verbose)
        return dag
    subprocess.run(cmd, cwd=job_dir, check=True)
    return dag
