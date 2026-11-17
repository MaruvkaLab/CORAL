"""The parallel stages must produce exactly what the serial ones produce: the
per-chromosome scan the same rows and triplet counts as the whole-genome scan,
and the flattened Fitch the same branch records as the recursive one."""
import csv
import gzip
import io
import random
from collections import defaultdict

from coral import parallel
from coral.multiple_species_mutation_extractor_manager import MultipleSpeciesMutationExtractor
from coral.multiple_species_utils import annotate_tree_with_indices

HEADER = ["chromosome", "position", "left", "right", "taxa0", "taxa1", "taxa2"]


def _pileup(chrom, n, seed, n_bams=2):
    """Random multi-species pileup lines; '.' means the sample matches the ref."""
    random.seed(seed)
    lines = []
    for i in range(n):
        ref = random.choice("ACGT")
        samples = "".join(f"\t5\t{random.choice(['.', '.', '.', 'T', 'G'])*5}\tIIIII"
                          for _ in range(n_bams))
        lines.append(f"{chrom}\t{100+i}\t{ref}{samples}")
    return lines


def _extractor(tmp_path, tree=None, mapping=None, species=None):
    return MultipleSpeciesMutationExtractor(
        pileup_file=None, output_dir=str(tmp_path), n_species=3, tree=tree,
        species_list=species or ["o", "a", "b"], mapping=mapping or {0: "o", 1: "a", 2: "b"},
        no_cache=True, verbose=False)


def test_scan_per_chromosome_matches_whole_genome(tmp_path):
    # two chromosomes: the serial scan sees them concatenated, the parallel scan
    # one stream each, so this also covers the chromosome boundary
    chroms = {"chr1": _pileup("chr1", 300, 1), "chr2": _pileup("chr2", 300, 2)}
    pileup = tmp_path / "in.pileup.gz"
    with gzip.open(pileup, "wt") as f:
        f.write("\n".join(chroms["chr1"] + chroms["chr2"]) + "\n")

    serial = _extractor(tmp_path)
    serial.pileup_file = str(pileup)
    csv_path, triplets_path = tmp_path / "m.csv.gz", tmp_path / "t.json"
    serial_triplets = serial._scan_pileup_serial(str(csv_path), str(triplets_path), HEADER)
    with gzip.open(csv_path, "rt") as f:
        serial_rows = list(csv.reader(f))[1:]

    parallel_rows, parallel_triplets = [], defaultdict(int)
    for chrom in ("chr1", "chr2"):
        stream = io.StringIO("\n".join(chroms[chrom]) + "\n")
        counts = parallel.scan_stream(_extractor(tmp_path), stream, parallel_rows.append)
        for context, n in counts.items():
            parallel_triplets[context] += n

    assert serial_rows == [[str(v) for v in row] for row in parallel_rows]
    assert serial_triplets == dict(parallel_triplets)
    assert serial_rows                                  # the fixture must detect something


def _fitch_csv(path, taxa, n_rows=300, seed=5):
    random.seed(seed)
    with gzip.open(path, "wt", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["chromosome", "position", "left", "right"] + taxa)
        for i in range(n_rows):
            writer.writerow(["chr1", 1000 + i, random.choice("ACGT"), random.choice("ACGT")]
                            + [random.choice("AACGT") for _ in taxa])


TREES = [
    "(((sp1,sp2),(sp3,sp4)),outg);",     # balanced
    "((((sp1,sp2),sp3),sp4),outg);",     # ladder
    "((sp1,sp2,sp3),sp4,outg);",         # multifurcating
    "(((sp2,sp1),(sp4,sp3)),outg);",     # same clades, children swapped
]


def test_flattened_fitch_matches_recursive(tmp_path):
    for i, newick in enumerate(TREES):
        tree, mapping, species = annotate_tree_with_indices(newick, "outg", verbose=False)
        taxa = [f"taxa{k}" for k in mapping if isinstance(k, int)]
        csv_path = tmp_path / f"fitch{i}.csv.gz"
        _fitch_csv(csv_path, taxa, seed=i)

        extractor = _extractor(tmp_path, tree=tree, mapping=mapping, species=species)
        extractor.n_species = len(species)
        recursive, recursive_ambiguous = extractor._fitch_serial(str(csv_path))
        flat, flat_ambiguous = parallel.run_fitch(str(csv_path), tree, mapping, fitch_jobs=1)

        assert set(recursive) == set(flat), newick
        for branch in recursive:                        # record order matters downstream
            assert list(recursive[branch]) == list(flat[branch]), newick
        assert recursive_ambiguous == flat_ambiguous, newick
        assert recursive                                # the fixture must find mutations


def test_job_count_never_exceeds_cores():
    assert parallel.cap_jobs(None, 8, 20) == 8          # auto
    assert parallel.cap_jobs(16, 8, 20) == 8            # explicit --*-jobs is clamped
    assert parallel.cap_jobs(2, 8, 20) == 2             # ...but may still reduce
    assert parallel.cap_jobs(None, 8, 3) == 3           # never more workers than tasks


def test_split_threads_allocates_every_core():
    assert parallel.split_threads(10, 4, 4) == [3, 3, 2, 2]     # remainder, not floor
    assert parallel.split_threads(16, 4, 4) == [4, 4, 4, 4]
    # when tasks queue, thread counts stay uniform so concurrent workers can't
    # oversubscribe: 2 workers x 8 threads live, not 4 x 8
    assert parallel.split_threads(16, 2, 4) == [8, 8, 8, 8]
