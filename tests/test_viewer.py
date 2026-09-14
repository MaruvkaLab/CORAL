"""coral view: inputs kept before cleanup, sister-mask tracing through fragment reads,
calls removed by the mask, region / size limits, and the written page."""
import csv
import gzip
import json
import random
import re
import subprocess
import sys
from types import SimpleNamespace

import numpy as np
import pysam
import pytest

from coral.viewer import (ViewerError, build_view_data, build_viewer, pack_sequence,
                          preserve_viewer_inputs, unpack_sequence)

REF, T1, T2 = "Ref_genome", "Taxon_one", "Taxon_two"
CONTIGS = [("chrA", 400), ("chrB", 120)]


def _reference_sequences():
    rng = random.Random(7)
    seqs = {name: "".join(rng.choice("ACGT") for _ in range(n)) for name, n in CONTIGS}
    seqs["chrB"] = seqs["chrB"][:50] + "NNNN" + seqs["chrB"][54:]       # 1-based 51..54
    return seqs


def _write_fasta(path, seqs):
    with open(path, "w") as fh:
        for name, seq in seqs.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")
    pysam.faidx(str(path))


def _write_bed(path, rows):
    with open(path, "w") as fh:
        for chrom, start0, end in rows:
            fh.write(f"{chrom}\t{start0}\t{end}\n")


def _write_calls(run, taxon, other, calls):
    mut = run / "Mutations"
    mut.mkdir(parents=True, exist_ok=True)
    with gzip.open(mut / f"{taxon}__{other}__{REF}__mutations.csv.gz", "wt", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["chromosome", "position", "mutation"])
        w.writerows(calls)


def _write_bam(path, reads):
    """reads: (name, contig, 0-based start, length, flag)."""
    header = {"HD": {"VN": "1.6", "SO": "unsorted"}, "SQ": [{"SN": n, "LN": l} for n, l in CONTIGS]}
    unsorted = str(path) + ".unsorted.bam"
    with pysam.AlignmentFile(unsorted, "wb", header=header) as out:
        for name, chrom, start0, length, flag in reads:
            a = pysam.AlignedSegment(out.header)
            a.query_name, a.flag, a.mapping_quality = name, flag, 60
            a.reference_id, a.reference_start = out.get_tid(chrom), start0
            a.cigarstring = f"{length}M"
            a.query_sequence = "A" * length
            a.query_qualities = pysam.qualitystring_to_array("I" * length)
            out.write(a)
    pysam.sort("-o", str(path), unsorted)
    pysam.index(str(path))


@pytest.fixture
def runs(tmp_path):
    """A masked run (with viewer inputs) and the same triplet run without masking."""
    run, cmp = tmp_path / f"{REF}__{T1}__{T2}_masked", tmp_path / f"{REF}__{T1}__{T2}_unmasked"
    seqs = _reference_sequences()
    genomes = tmp_path / "genomes"
    genomes.mkdir()
    _write_fasta(genomes / "ref.fasta", seqs)
    _write_bed(genomes / "ref.bed", [("chrA", 150, 160)])                 # 1-based 151..160
    _write_bed(genomes / "t1.bed", [("sA", 20, 21), ("sA", 338, 339)])    # 1-based 21 and 339
    _write_bed(genomes / "t2.bed", [])

    preserve_viewer_inputs(
        str(run),
        SimpleNamespace(name=REF, accession="GCF_REF", fasta_path=str(genomes / "ref.fasta")),
        [SimpleNamespace(name=T1, accession="GCF_T1"), SimpleNamespace(name=T2, accession="GCF_T2")],
        {REF: str(genomes / "ref.bed"), T1: str(genomes / "t1.bed"), T2: str(genomes / "t2.bed")},
        parameters={"repeat_mask": True, "repeat_masker": "windowmasker", "mapq": 60, "low_mapq": 1,
                    "aligner_name": "bwa", "nested": {"dropped": True}},
        verbose=False)

    _write_calls(run, T1, T2, [("chrA", 120, "A[C>T]G"), ("chrA", 130, "C[A>G]T"),
                               ("chrA", 210, "G[T>C]A"), ("chrA", 300, "T[G>A]T"), ("chrB", 10, "A[C>A]A")])
    _write_calls(run, T2, T1, [("chrA", 50, "A[A>G]A")])
    (run / "BAMs").mkdir()
    _write_bam(run / "BAMs" / f"{T1}_to_{REF}.bam", [
        ("sA_1_50", "chrA", 99, 50, 0),        # 100..149 fwd: call 120 <- sister 21 (masked), 130 <- 31
        ("sA_200_249", "chrA", 124, 50, 0),    # 125..174 fwd: call 130 <- sister 205
        ("sA_300_349", "chrA", 199, 50, 16),   # 200..249 rev: call 210 <- sister 349-10 = 339 (masked)
        ("sA_1_50", "chrA", 279, 50, 256),     # secondary over call 300 from masked sister 21: ignored
    ])
    _write_bam(run / "BAMs" / f"{T2}_to_{REF}.bam", [])

    _write_calls(cmp, T1, T2, [("chrA", 120, "A[C>T]G"), ("chrA", 130, "C[A>G]T"), ("chrA", 155, "C[C>T]A"),
                               ("chrA", 161, "G[A>C]C"), ("chrA", 210, "G[T>C]A"), ("chrA", 300, "T[G>A]T"),
                               ("chrA", 380, "T[C>G]A"), ("chrB", 10, "A[C>A]A")])
    _write_calls(cmp, T2, T1, [("chrA", 50, "A[A>G]A"), ("chrA", 90, "A[T>A]A")])
    return SimpleNamespace(run=run, cmp=cmp, seqs=seqs)


def _positions(packed):
    return np.cumsum(packed["d"]).tolist()


def _mask(contig):
    d, out, p = contig["mask"], [], 0
    for i in range(0, len(d), 2):
        p += d[i]
        out.append((p, p + d[i + 1] - 1))
    return out


def test_pack_sequence_round_trip_with_non_acgt():
    seq = "ACGTNNacgtRYA"
    packed, runs = pack_sequence(seq)
    assert runs == [4, 2, 10, 2]
    assert unpack_sequence(packed, len(seq), runs) == "ACGTNNACGTNNA"


def test_preserve_keeps_indexed_fasta_masks_and_manifest(runs):
    inputs = runs.run / "Viewer" / "inputs"
    manifest = json.loads((inputs / "manifest.json").read_text())
    assert manifest["taxa"] == [T1, T2] and manifest["reference"] == REF
    assert manifest["masks"] == {REF: f"mask.{REF}.bed", T1: f"mask.{T1}.bed", T2: f"mask.{T2}.bed"}
    assert "nested" not in manifest["parameters"]
    with pysam.FastaFile(str(inputs / manifest["fasta"])) as fa:
        assert fa.fetch("chrA") == runs.seqs["chrA"]


def test_view_data_calls_flags_and_removed(runs):
    data = build_view_data(str(runs.run), compare=str(runs.cmp))
    assert data["branches"] == [T1, T2]
    meta = data["meta"]
    assert meta["has_mask"] and meta["has_sister_masks"] and meta["compare"] == runs.cmp.name
    chrA, chrB = data["contigs"]
    assert (chrA["name"], chrA["start"], chrA["end"]) == ("chrA", 1, 400)

    assert unpack_sequence(chrA["seq"], 400) == runs.seqs["chrA"]
    n_offsets = [chrB["n"][0] - chrB["start"], chrB["n"][1]]
    assert chrB["n"] == [51, 4]
    assert unpack_sequence(chrB["seq"], 120, n_offsets) == runs.seqs["chrB"]
    assert _mask(chrA) == [(151, 160)] and chrB["mask"] == []

    kept = chrA["kept"][T1]
    assert _positions(kept) == [120, 130, 210, 300]
    assert kept["f"] == [1, 0, 1, 0]           # 120 and 210 only from masked sister bases; 300 only a secondary
    assert [data["contexts"][i] for i in kept["c"]] == ["A[C>T]G", "C[A>G]T", "G[T>C]A", "T[G>A]T"]
    assert chrB["kept"][T1]["f"] == [0]          # no reads -> not flagged

    removed = chrA["removed"][T1]
    assert _positions(removed) == [155, 161, 380]
    assert removed["f"] == [0, 0, 1]             # in / next to the reference mask; elsewhere -> sister mask
    assert _positions(chrA["removed"][T2]) == [90] and chrA["removed"][T2]["f"] == [1]
    assert _positions(chrA["kept"][T2]) == [50] and chrA["kept"][T2]["f"] == [0]


def test_region_limits_sequence_mask_and_calls(runs):
    data = build_view_data(str(runs.run), compare=str(runs.cmp), region="chrA:100-250")
    (chrA,) = data["contigs"]
    assert (chrA["start"], chrA["end"], chrA["length"]) == (100, 250, 400)
    assert unpack_sequence(chrA["seq"], 151) == runs.seqs["chrA"][99:250]
    assert _positions(chrA["kept"][T1]) == [120, 130, 210]
    assert _positions(chrA["removed"][T1]) == [155, 161]
    assert _mask(chrA) == [(151, 160)]


def test_limits_and_errors(runs, tmp_path):
    with pytest.raises(ViewerError, match="--region"):
        build_view_data(str(runs.run), max_genome_bp=100)
    with pytest.raises(ViewerError, match="not found in the reference"):
        build_view_data(str(runs.run), region="chrZ:1-10")
    with pytest.raises(ViewerError, match="--viewer"):
        build_view_data(str(tmp_path / "no_such_run"))

    manifest_path = runs.run / "Viewer" / "inputs" / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["parameters"]["repeat_mask"] = False
    manifest_path.write_text(json.dumps(manifest))
    with pytest.raises(ViewerError, match="without --repeat-mask"):
        build_view_data(str(runs.run), compare=str(runs.cmp))


def test_build_viewer_embeds_data_in_page(runs):
    out = build_viewer(str(runs.run), compare=str(runs.cmp), verbose=False)
    assert out.endswith(f"Viewer/{runs.run.name}.html")
    html = open(out, encoding="utf-8").read()
    assert "__DATA__" not in html
    payload = re.search(r'<script id="data" type="application/json">(.*?)</script>', html, re.S).group(1)
    assert json.loads(payload.replace("<\\/", "</")) == build_view_data(str(runs.run), compare=str(runs.cmp))


def test_cli_view(runs, tmp_path):
    out = tmp_path / "page.html"
    ok = subprocess.run([sys.executable, "-m", "coral", "view", str(runs.run), "--compare", str(runs.cmp),
                         "--region", "chrA", "--output", str(out), "--quiet"], capture_output=True, text=True)
    assert ok.returncode == 0, ok.stderr
    assert out.exists()
    bad = subprocess.run([sys.executable, "-m", "coral", "view", str(tmp_path / "missing")], capture_output=True, text=True)
    assert bad.returncode == 1 and "--viewer" in bad.stderr
