import textwrap

import numpy as np

from coral.repeat_masker import RepeatMask, _merge


def _merged(ivs):
    s, e = _merge(np.asarray(ivs, dtype=np.int64).reshape(-1, 2))
    return list(zip(s.tolist(), e.tolist()))


def test_merge_touching_and_overlapping():
    assert _merged([(10, 20), (21, 30)]) == [(10, 30)]      # adjacent -> merged
    assert _merged([(10, 20), (15, 25)]) == [(10, 25)]      # overlapping
    assert _merged([(50, 60), (10, 20)]) == [(10, 20), (50, 60)]  # sorted, disjoint


def test_contains_and_fraction():
    m = RepeatMask({"chr1": [(100, 200), (300, 350)]})
    assert m
    assert not m.contains("chr1", 99)
    assert m.contains("chr1", 100) and m.contains("chr1", 200)   # inclusive ends
    assert not m.contains("chr1", 201)
    assert m.contains("chr1", 325)
    assert not m.contains("chrX", 150)                            # unknown contig
    # span [150,260]: masked [150,200] = 51 bases of 111
    assert m.masked_bases("chr1", 150, 260) == 51
    assert abs(m.masked_fraction("chr1", 150, 260) - 51 / 111) < 1e-9


def test_for_contig_slices_to_one_contig():
    m = RepeatMask({"chr1": [(100, 200)], "chr2": [(10, 90)]})
    sub = m.for_contig("chr1")
    assert sub.contains("chr1", 150)
    assert not sub.contains("chr2", 50)          # chr2 not in the slice
    assert m.for_contig("chrZ").n_intervals == 0  # unknown contig -> empty mask


def test_empty_mask_passes_everything():
    m = RepeatMask()
    assert not m
    assert not m.contains("chr1", 10)
    assert m.masked_fraction("chr1", 1, 150) == 0.0


def test_from_out_parses_repeatmasker(tmp_path):
    out = tmp_path / "ref.fasta.out"
    out.write_text(textwrap.dedent("""\
       SW   perc perc perc  query   position in query   matching  repeat        position in repeat
    score   div. del. ins.  sequence begin end (left)   repeat    class/family  begin end (left) ID

      463   1.3  0.6  1.7  chr1       105   220 (900) +  (TAACCC)n Simple_repeat     1  116  (0)  1
     1200   8.4  0.0  0.0  chr1       500   640 (480) C  L1MEg     LINE/L1        (0) 200   61   2
      300   2.0  0.0  0.0  chr2        10    90 (700) +  (CA)n     Simple_repeat     1   80  (0)  3
    """))
    m = RepeatMask.from_out(str(out))
    assert m.contains("chr1", 105) and m.contains("chr1", 220)
    assert m.contains("chr1", 500) and m.contains("chr1", 640)
    assert not m.contains("chr1", 300)          # gap between the two chr1 repeats
    assert m.contains("chr2", 50)
    assert m.n_intervals == 3


def test_from_windowmasker_shifts_to_1based(tmp_path):
    iv = tmp_path / "ref.fasta.wm_intervals"
    # WindowMasker interval coords are 0-based inclusive -> our 1-based +1
    iv.write_text(">chr1\n922 - 928\n1000 - 1010\n>chr2\n5 - 9\n")
    m = RepeatMask.from_windowmasker(str(iv))
    assert m.contains("chr1", 923) and m.contains("chr1", 929)   # 922..928 -> 923..929
    assert not m.contains("chr1", 930)
    assert m.contains("chr2", 6) and m.contains("chr2", 10)
    assert m.n_intervals == 3


def test_bed_roundtrip(tmp_path):
    m = RepeatMask({"chr1": [(100, 200)], "chr2": [(10, 90)]})
    bed = tmp_path / "r.bed"
    m.to_bed(str(bed))
    # BED is 0-based half-open: [99,200) and [9,90)
    assert bed.read_text().splitlines()[0] == "chr1\t99\t200"
    back = RepeatMask.from_bed(str(bed))
    assert back.contains("chr1", 100) and back.contains("chr1", 200)
    assert not back.contains("chr1", 201)
    assert back.contains("chr2", 10) and back.contains("chr2", 90)
