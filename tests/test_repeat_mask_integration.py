"""The two behavioural changes repeat-masking makes: drop calls at reference repeat
positions (scan_pileup) and don't generate pseudo-reads from species repeats
(generate_fragment_fastq)."""
from coral.mutation_extractor_manager import scan_pileup
from coral.repeat_masker import RepeatMask
from coral.genome_manager import Genome


# raw 3-way mpileup cols: chrom pos ref depth1 bases1 quals1 depth2 bases2 quals2
def _pline(pos, ref, b1, b2):
    return f"chr1\t{pos}\t{ref}\t1\t{b1}\tI\t1\t{b2}\tI\n"

# center (pos 101) is a taxon-2 substitution C>T; flanks match in both taxa
PILEUP = [
    _pline(100, "A", ".", "."),
    _pline(101, "C", ".", "T"),
    _pline(102, "G", ".", "."),
]


def test_unmasked_calls_the_mutation():
    mut1, mut2, _, _ = scan_pileup(iter(PILEUP))
    assert mut2 == {"A[C>T]G": 1}
    assert mut1 == {}


def test_masking_the_center_drops_the_call():
    mask = RepeatMask({"chr1": [(101, 101)]})
    mut1, mut2, _, _ = scan_pileup(iter(PILEUP), ref_mask=mask)
    assert mut2 == {}                       # position removed -> no call there


def test_masking_a_flank_also_drops_the_call():
    # removing a flank breaks the 3-position window, so the center can't be called
    mask = RepeatMask({"chr1": [(100, 100)]})
    _, mut2, _, _ = scan_pileup(iter(PILEUP), ref_mask=mask)
    assert mut2 == {}


def test_mask_elsewhere_keeps_the_call():
    mask = RepeatMask({"chr1": [(500, 600)], "chr2": [(1, 100)]})
    _, mut2, _, _ = scan_pileup(iter(PILEUP), ref_mask=mask)
    assert mut2 == {"A[C>T]G": 1}


def test_fragment_generation_skips_repeat_reads(tmp_path):
    fa = tmp_path / "sp.fasta"
    fa.write_text(">chr1\n" + "ACGT" * 75 + "\n")     # 300 bp
    g = Genome("sp", "ACC", str(tmp_path), fasta_path=str(fa), verbose=False)

    # length 150, offset 75 -> fragments (1-150), (76-225), (151-300)
    mask = RepeatMask({"chr1": [(1, 170)]})           # first ~170 bp is repeat
    out = tmp_path / "sp.fastq"
    g.generate_fragment_fastq(length=150, offset=75, output_fastq=str(out),
                              force=True, repeat_mask=mask, mask_frac=0.5)
    names = [ln[1:].strip() for ln in out.read_text().splitlines() if ln.startswith("@")]
    # (1-150) fully masked and (76-225) 63% masked are dropped; (151-300) kept
    assert names == ["chr1_151_300"]


def test_fragment_generation_unmasked_keeps_all(tmp_path):
    fa = tmp_path / "sp.fasta"
    fa.write_text(">chr1\n" + "ACGT" * 75 + "\n")
    g = Genome("sp", "ACC", str(tmp_path), fasta_path=str(fa), verbose=False)
    out = tmp_path / "sp.fastq"
    g.generate_fragment_fastq(length=150, offset=75, output_fastq=str(out), force=True)
    names = [ln[1:].strip() for ln in out.read_text().splitlines() if ln.startswith("@")]
    assert names == ["chr1_1_150", "chr1_76_225", "chr1_151_300"]
