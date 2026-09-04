"""Site classification in the multi-species scan, and reading back the tree
PHYLIP writes. taxa0 is the outgroup: it has no BAM, so its base is the
reference base and the CSV row carries it as the fifth fixed column."""
from coral.multiple_species_mutation_extractor_manager import MultipleSpeciesMutationExtractor
from coral.multiple_species_utils import (annotate_tree_with_indices, parse_phylip_edges,
                                          phylip_interior_clades, read_phylip_outtree_newick,
                                          tree_from_phylip_outtree)

# 3-line window; ref at the centre is C, flanks A and G match in both samples
def _window(extractor, sample1, sample2):
    lines = [f"chr1\t100\tA\t5\t.....\tIIIII\t5\t.....\tIIIII",
             f"chr1\t101\tC\t5\t{sample1*5}\tIIIII\t5\t{sample2*5}\tIIIII",
             f"chr1\t102\tG\t5\t.....\tIIIII\t5\t.....\tIIIII"]
    return [extractor._parse_line(line) for line in lines]


def _extractor(tmp_path):
    return MultipleSpeciesMutationExtractor(
        pileup_file=None, output_dir=str(tmp_path), n_species=3, tree=None,
        species_list=["o", "a", "b"], mapping={0: "o", 1: "a", 2: "b"},
        no_cache=True, verbose=False)


def test_ingroup_varies_gives_a_row_and_an_opportunity(tmp_path):
    e = _extractor(tmp_path)
    row, triplet = e._detect_site(_window(e, ".", "T"))
    assert row == ["chr1", "101", "A", "G", "C", "C", "T"]
    assert triplet == "ACG"
    assert len(row) == len(["chromosome", "position", "left", "right", "taxa0", "taxa1", "taxa2"])


def test_invariant_site_is_an_opportunity_but_not_a_row(tmp_path):
    e = _extractor(tmp_path)
    row, triplet = e._detect_site(_window(e, ".", "."))
    assert row is None
    assert triplet == "ACG"


def test_site_with_no_sample_at_the_reference_base_is_dropped(tmp_path):
    # a single outgroup cannot polarize these, so they count for neither the
    # numerator nor the denominator
    e = _extractor(tmp_path)
    assert e._detect_site(_window(e, "T", "T")) == (None, None)   # uniform, differs from ref
    assert e._detect_site(_window(e, "T", "G")) == (None, None)   # tri-allelic


def test_unconserved_flank_drops_the_site(tmp_path):
    e = _extractor(tmp_path)
    buffer = _window(e, ".", "T")
    buffer[0] = e._parse_line("chr1\t100\tA\t5\t.....\tIIIII\t5\tGGGGG\tIIIII")
    assert e._detect_site(buffer) == (None, None)


# dnapars writes the newick alone; with a user tree ('U') it prefixes a count line
OUTTREE = "((taxa1:0.1,taxa2:0.1):0.05,(taxa3:0.1,\ntaxa4:0.1):0.05,taxa0:0.1);\n"
OUTFILE = """
between        and            length
-------        ---            ------
   6          taxa0           0.10000
   6             7            0.05000
   7          taxa1           0.10000
   7          taxa2           0.10000
   6             8            0.05000
   8          taxa3           0.10000
   8          taxa4           0.10000

remember: this is an unrooted tree!
"""


def _phylip_run(tmp_path, name, outtree):
    (tmp_path / f"{name}.outtree").write_text(outtree)
    (tmp_path / f"{name}.outfile").write_text(OUTFILE)
    return str(tmp_path / f"{name}.outtree")


def test_outtree_reads_with_and_without_a_leading_count_line(tmp_path):
    plain = read_phylip_outtree_newick(_phylip_run(tmp_path, "plain", OUTTREE), verbose=False)
    counted = read_phylip_outtree_newick(_phylip_run(tmp_path, "counted", "1\n" + OUTTREE),
                                         verbose=False)
    assert plain == counted                             # the count line must not survive
    assert plain.startswith("(") and plain.endswith(";")
    assert "\n" not in plain                            # PHYLIP wraps at ~70 columns


def test_phylip_tree_keeps_the_csv_column_indices(tmp_path):
    _, mapping, _ = annotate_tree_with_indices("(((sp1,sp2),(sp3,sp4)),outg);", "outg",
                                               verbose=False)
    tree = tree_from_phylip_outtree(_phylip_run(tmp_path, "run", OUTTREE), mapping, "outg")

    # leaf indices must come from the mapping, not a re-sort, or they desync
    # from the taxaK columns already written to matching_bases.csv.gz
    assert {leaf.name: leaf.index for leaf in tree.iter_leaves()} == \
        {name: idx for name, idx in mapping.items() if isinstance(name, str)}
    assert tree.get_leaves_by_name("outg")               # rooted on the outgroup


def test_interior_nodes_are_named_from_the_outfile(tmp_path):
    _, mapping, _ = annotate_tree_with_indices("(((sp1,sp2),(sp3,sp4)),outg);", "outg",
                                               verbose=False)
    tree = tree_from_phylip_outtree(_phylip_run(tmp_path, "run", OUTTREE), mapping, "outg")
    interior = {node.custom_name for node in tree.traverse() if not node.is_leaf()}
    assert {"6", "7", "8"} <= interior                   # PHYLIP's own numbering

    clades = phylip_interior_clades(parse_phylip_edges(str(tmp_path / "run.outfile")), "taxa0")
    assert clades["7"] == frozenset({"taxa1", "taxa2"})
    assert clades["8"] == frozenset({"taxa3", "taxa4"})
