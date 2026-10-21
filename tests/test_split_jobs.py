"""run_multi --split-jobs: off by default, and when on, the per-species stages
must make the same calls the in-process pipeline makes."""
import os
import subprocess
import sys

import pytest

from coral import cli, split_jobs
from coral import pipeline as pipeline_module
from coral.pipeline import MultiSpeciesMutationPipeline

NEWICK = "(Out|GCF_0,(A|GCF_1,(B|GCF_2,C_c|GCF_3)));"


class FakeGenome:
    calls = []

    def __init__(self, name, accession, output_dir, fasta_path=None, no_cache=False, verbose=True):
        FakeGenome.calls.append(("init", name, accession, output_dir, no_cache))
        self.name = name
        self.accession = accession
        self.output_dir = os.path.join(output_dir, name)
        self.fasta_path = os.path.join(self.output_dir, f"{name}.fasta")
        self.fastq_path = None

    def download(self):
        FakeGenome.calls.append(("download", self.name))
        os.makedirs(self.output_dir, exist_ok=True)
        if not os.path.exists(self.fasta_path):
            with open(self.fasta_path, "w") as f:
                f.write(">c\nACGT\n")

    def index(self, aligner="bwa"):
        FakeGenome.calls.append(("index", self.name, aligner))
        for ext in split_jobs.INDEX_EXTS.get(aligner, []) + ["fai"]:
            open(f"{self.fasta_path}.{ext}", "w").close()

    def generate_fragment_fastq(self, length=150, output_fastq=None, offset=75, force=False,
                                repeat_mask=None, mask_frac=0.5):
        out = output_fastq or os.path.join(self.output_dir, f"{self.name}.fastq")
        if not os.path.exists(out) or force:
            FakeGenome.calls.append(("fragment", self.name, length, offset, repeat_mask))
            with open(out, "w") as f:
                f.write("@r1\nACGT\n+\nIIII\n")
        self.fastq_path = out
        return out


class FakeAligner:
    calls = []

    def __init__(self, species_genome, reference_genome, base_output_dir, aligner_cmd=None,
                 aligner_name=None, no_cache=False, cores=None, verbose=True):
        FakeAligner.calls.append(("init", species_genome.name, reference_genome.name,
                                  reference_genome.fasta_path, species_genome.fastq_path,
                                  base_output_dir, aligner_cmd, aligner_name, no_cache))
        self.cores = cores
        self.species = species_genome.name
        self.fastq = species_genome.fastq_path
        bam_dir = os.path.join(base_output_dir, "BAMs")
        os.makedirs(bam_dir, exist_ok=True)
        self.raw_bam = os.path.join(bam_dir, f"{self.species}_raw.bam")
        self.final_bam = os.path.join(bam_dir, f"{self.species}.bam")
        self.filter_stats = None

    def align_disk_cached(self, **kwargs):
        FakeAligner.calls.append(("disk_cached", self.species, sorted(kwargs.items())))
        open(self.raw_bam, "w").close()
        open(self.final_bam, "w").close()

    def align_streamed(self, **kwargs):
        FakeAligner.calls.append(("streamed", self.species, sorted(kwargs.items())))
        open(self.final_bam, "w").close()


@pytest.fixture
def fakes(monkeypatch):
    FakeGenome.calls = []
    FakeAligner.calls = []
    for module in (pipeline_module, split_jobs):
        monkeypatch.setattr(module, "Genome", FakeGenome)
        monkeypatch.setattr(module, "Aligner", FakeAligner)
    monkeypatch.setattr(split_jobs, "_check_raw_bam", lambda aligner, threads: None)
    monkeypatch.setattr(split_jobs, "_check_final_bam", lambda aligner: None)


def _pipeline(out, **kw):
    return MultiSpeciesMutationPipeline(newick_tree=NEWICK, base_output_dir=str(out),
                                        run_id="r", aligner_name="bwa", cores=2,
                                        no_parallel=True, verbose=False, **kw)


@pytest.mark.parametrize("streamed", [False, True])
def test_split_stages_make_the_in_process_calls(tmp_path, fakes, streamed):
    whole = _pipeline(tmp_path / "whole", streamed=streamed, mapq=50, low_mapq=2)
    whole.parse_and_annotate_tree()
    whole.download_index_and_fragment()
    whole.align_species_to_outgroup()
    genome_calls, aligner_calls = FakeGenome.calls, FakeAligner.calls

    FakeGenome.calls, FakeAligner.calls = [], []
    split_out = tmp_path / "split"
    split_jobs.run_prepare_stage(_pipeline(split_out, streamed=streamed, mapq=50, low_mapq=2))
    for i in range(3):
        split_jobs.run_align_stage(_pipeline(split_out, streamed=streamed, mapq=50, low_mapq=2),
                                   i, threads=3)

    def norm(calls):
        return [tuple(str(x).replace("/whole", "/split").replace(".fastq.tmp", ".fastq")
                      for x in c) for c in calls]

    # Every genome is downloaded and prepared once, with the same arguments;
    # the split run only adds the outgroup lookups the align jobs do.
    split_genome = [c for c in FakeGenome.calls if not (c[0] == "init" and c[1] == "Out")]
    whole_genome = [c for c in genome_calls if not (c[0] == "init" and c[1] == "Out")]
    assert sorted(norm(split_genome)) == sorted(norm(whole_genome))
    assert norm(FakeAligner.calls) == norm(aligner_calls)
    for species in ["A", "B", "C_c"]:
        assert (split_out / "r" / species / f"{species}.fastq").exists()
        assert not (split_out / "r" / species / f"{species}.fastq.tmp").exists()


def test_align_stage_refuses_an_unprepared_outgroup(tmp_path, fakes):
    with pytest.raises(FileNotFoundError, match="prepare"):
        split_jobs.run_align_stage(_pipeline(tmp_path), 0, threads=1)


def test_align_stage_index_out_of_range(tmp_path, fakes):
    split_jobs.run_prepare_stage(_pipeline(tmp_path))
    with pytest.raises(IndexError):
        split_jobs.run_align_stage(_pipeline(tmp_path), 3, threads=1)


def test_threads_leave_one_cpu():
    assert split_jobs.align_threads(4) == 3
    assert split_jobs.align_threads(1) == 1
    with pytest.raises(ValueError):
        split_jobs.align_threads(0)


def test_truncated_raw_bam_is_set_aside(tmp_path, monkeypatch):
    fastq = tmp_path / "x.fastq"
    fastq.write_text("@a\nAC\n+\nII\n@b\nAC\n+\nII\n")
    raw, final = tmp_path / "x_raw.bam", tmp_path / "x.bam"
    raw.write_text("raw")
    final.write_text("final")

    class A:
        pass
    aligner = A()
    aligner.raw_bam, aligner.final_bam, aligner.fastq = str(raw), str(final), str(fastq)

    monkeypatch.setattr(split_jobs.subprocess, "run",
                        lambda *a, **k: subprocess.CompletedProcess(a, 0, stdout="2\n"))
    split_jobs._check_raw_bam(aligner, 1)
    assert raw.exists()

    monkeypatch.setattr(split_jobs.subprocess, "run",
                        lambda *a, **k: subprocess.CompletedProcess(a, 0, stdout="1\n"))
    with pytest.raises(RuntimeError, match="did not finish"):
        split_jobs._check_raw_bam(aligner, 1)
    assert not raw.exists() and (tmp_path / "x_raw.bam.incomplete").exists()
    assert not final.exists() and (tmp_path / "x.bam.incomplete").exists()


def _run_cli(monkeypatch, argv):
    monkeypatch.setattr(sys, "argv", ["coral"] + argv)
    cli.main()


def test_without_the_flag_run_multi_is_unchanged(tmp_path, monkeypatch):
    ran = []
    monkeypatch.setattr(MultiSpeciesMutationPipeline, "run", lambda self: ran.append(self))

    def no_split(*a, **k):
        raise AssertionError("split_jobs must not run without --split-jobs")
    monkeypatch.setattr(split_jobs, "dispatch", no_split)
    _run_cli(monkeypatch, ["run_multi", "--newick-tree", NEWICK, "--output", str(tmp_path),
                           "--cores", "2"])
    assert len(ran) == 1
    assert not (tmp_path / "multi_species_run" / "split_jobs").exists()


def test_main_stage_is_plain_run(tmp_path, monkeypatch):
    ran = []
    monkeypatch.setattr(MultiSpeciesMutationPipeline, "run", lambda self: ran.append(self.cores))
    _run_cli(monkeypatch, ["run_multi", "--newick-tree", NEWICK, "--output", str(tmp_path),
                           "--split-jobs", "--split-stage", "main", "--cores", "3"])
    assert ran == [3]


@pytest.fixture
def submitted(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(split_jobs, "check_phylip_available", lambda cmd: True)
    argv = ["run_multi", "--newick-tree", NEWICK, "--output", "out", "--run-id", "r",
            "--aligner-name", "bwa", "--mapq", "60", "--split-jobs", "--no-submit",
            "--job-cores", "4", "--job-extra", "request_memory = 12000",
            "--main-job-extra", "+Group = \"x\""]
    _run_cli(monkeypatch, argv)
    return tmp_path / "out" / "r" / "split_jobs", argv


def test_submit_writes_dag_and_jobs(submitted):
    job_dir, argv = submitted
    dag = (job_dir / "r.dag").read_text()
    assert "JOB prepare" in dag and "JOB main" in dag
    assert [line.split()[1] for line in dag.splitlines() if line.startswith("JOB align_")] \
        == ["align_0000", "align_0001", "align_0002"]
    assert "# align_0002 = C_c" in dag
    assert "MAXJOBS align 20" in dag

    for name in ["prepare", "align", "main"]:
        sub = (job_dir / f"{name}.sub").read_text()
        assert "transfer" not in sub.lower()
        lines = sub.splitlines()
        assert lines[-1] == "queue"
    align = (job_dir / "align.sub").read_text().splitlines()
    assert "request_cpus   = 4" in align
    assert align[-2] == "request_memory = 12000"   # extra lines come last and win
    main = (job_dir / "main.sub").read_text().splitlines()
    assert main[-2] == '+Group = "x"'
    assert "request_cpus   = 8" in main            # 2 x the --cores 4 passed to it

    script = job_dir / "run_stage.sh"
    subprocess.run(["bash", "-n", str(script)], check=True)
    assert os.access(script, os.X_OK)


def test_stage_script_passes_argv_verbatim(tmp_path):
    fake_python = tmp_path / "py"
    fake_python.write_text('#!/bin/bash\nprintf "%s\\n" "$@"\n')
    fake_python.chmod(0o755)
    argv = ["run_multi", "--newick-tree", "(a b|1,'c'|2);", "--output", "dir with space"]
    script = tmp_path / "run_stage.sh"
    script.write_text(split_jobs.stage_script(argv, str(tmp_path), str(fake_python),
                                              [("FOO", "a b")], ["--cores", "4"]))
    out = subprocess.run(["bash", str(script), "align", "7"], capture_output=True,
                         text=True, check=True).stdout.splitlines()
    assert out == ["-m", "coral"] + argv + ["--split-stage", "align", "--split-index", "7"]
    out = subprocess.run(["bash", str(script), "main"], capture_output=True,
                         text=True, check=True).stdout.splitlines()
    assert out[-4:] == ["--split-stage", "main", "--cores", "4"]


def test_submit_refuses_no_cache_and_a_live_lock(tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(split_jobs, "check_phylip_available", lambda cmd: True)
    base = ["run_multi", "--newick-tree", NEWICK, "--output", "out", "--run-id", "r",
            "--split-jobs", "--no-submit"]
    with pytest.raises(SystemExit):
        _run_cli(monkeypatch, base + ["--no-cache"])
    assert "--no-cache" in capsys.readouterr().err

    job_dir = tmp_path / "out" / "r" / "split_jobs"
    job_dir.mkdir(parents=True)
    (job_dir / "r.dag.lock").write_text("")
    with pytest.raises(SystemExit):
        _run_cli(monkeypatch, base)
    assert "lock" in capsys.readouterr().err
