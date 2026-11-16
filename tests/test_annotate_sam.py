"""annotate_sam must let the production continuity filter be reproduced post hoc:
for any MAPQ threshold, a read is kept by with_continuity_filter_sam exactly when
its MAPQ >= threshold or its ZN tag >= 0."""
import subprocess
from io import TextIOWrapper

import pysam
import pytest

from coral.alignment_manager import annotate_sam, with_continuity_filter_sam

LEN = 150

# Records in aligner output order (= source-genome fragment order), not sorted.
# (name, target contig, 1-based target start, MAPQ, flag)
RECORDS = [
    ("f0", "chr1", 1000, 60, 0),
    ("f1", "chr1", 1075, 20, 0),     # overlaps f0
    ("f2", "chr1", 5000, 30, 0),     # disjoint from f1; next surviving neighbour is f4
    ("f3", "chr1", 8000, 0, 0),      # below low_mapq: dropped before neighbours are set
    ("f4", "chr1", 5075, 10, 0),     # overlaps f2
    ("f5", "chr1", 9000, 5, 0),      # isolated
    ("f5", "chr1", 5150, 60, 2048),  # supplementary piece of f5: same group, overlaps neighbour f4
    ("f6", None, None, 0, 4),        # unmapped
    ("f7", "chrUn_x", 9000, 60, 0),  # alt contig: still a neighbour, but never written
    ("f8", "chr1", 9100, 3, 0),      # previous neighbour is on another contig
]

# Keyed (name, is_supplementary). f4 is rescued by f5's supplementary (60), which in
# turn overlaps f4 (10); f5's primary overlaps nothing.
EXPECTED_ZN = {("f0", False): 20, ("f1", False): 60, ("f2", False): 10, ("f4", False): 60,
               ("f5", False): -1, ("f5", True): 10, ("f8", False): -1}


def _write_bam(path):
    header = {"HD": {"VN": "1.6", "SO": "unsorted"},
              "SQ": [{"SN": "chr1", "LN": 20000}, {"SN": "chrUn_x", "LN": 20000}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for name, contig, start, mapq, flag in RECORDS:
            read = pysam.AlignedSegment(bam.header)
            read.query_name = name
            read.flag = flag
            read.query_sequence = "A" * LEN
            read.query_qualities = pysam.qualitystring_to_array("I" * LEN)
            if not flag & 4:
                read.reference_name = contig
                read.reference_start = start - 1
                read.cigarstring = f"{LEN}M"
                read.mapping_quality = mapq
            bam.write(read)


def _read_sam(path):
    with pysam.AlignmentFile(str(path), "r") as sam:
        return {(r.query_name, r.is_supplementary): r for r in sam}


@pytest.fixture
def bam(tmp_path):
    path = tmp_path / "raw.bam"
    _write_bam(path)
    return path


@pytest.fixture
def annotated(bam, tmp_path):
    out = tmp_path / "annotated.sam"
    annotate_sam(str(bam), str(out), verbose=False)
    return _read_sam(out)


def test_writes_mapped_main_contig_reads_with_expected_tags(annotated):
    assert {name: r.get_tag("ZN") for name, r in annotated.items()} == EXPECTED_ZN


@pytest.mark.parametrize("threshold", [1, 3, 5, 10, 20, 30, 60, 61])
def test_reproduces_continuity_filter(bam, annotated, tmp_path, threshold):
    out = tmp_path / f"filtered_{threshold}.sam"
    # Fed exactly as align_disk_cached does: a samtools view pipe. The production
    # filter's fetch() lacks until_eof, so an unindexed BAM path would raise.
    proc = subprocess.Popen(["samtools", "view", "-h", str(bam)], stdout=subprocess.PIPE)
    with_continuity_filter_sam(TextIOWrapper(proc.stdout), str(out),
                               mapq_threshold=threshold, verbose=False)
    proc.wait()
    derived = {name for name, r in annotated.items()
               if r.mapping_quality >= threshold or r.get_tag("ZN") >= 0}
    assert set(_read_sam(out)) == derived


@pytest.mark.parametrize("continuity", [False, True])
def test_align_disk_cached_keeps_every_filtered_read(bam, tmp_path, continuity):
    """Regression: filter_sam wrote through an unflushed text buffer and
    align_disk_cached closed the pipe under it, silently dropping the last reads.
    With this small input every read would still be in that buffer."""
    from types import SimpleNamespace
    from coral.alignment_manager import Aligner

    aligner = Aligner(SimpleNamespace(name="sp", fasta_path="sp.fa", fastq_path="sp.fq"),
                      SimpleNamespace(name="ref", fasta_path="ref.fa"),
                      str(tmp_path / "run"), cores=1, verbose=False)
    aligner.raw_bam = str(bam)
    aligner.align_disk_cached(mapq=10, low_mapq=1, continuity=continuity)
    with pysam.AlignmentFile(aligner.final_bam) as out:
        written = sum(1 for _ in out)
    assert written == aligner.filter_stats["kept_reads"] > 0
