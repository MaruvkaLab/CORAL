"""site_intervals must agree with production quality_check at every threshold."""
import os
import random
from collections import defaultdict

import pytest

from coral.annotation import ALWAYS_KEPT, site_intervals, split_reads
from coral.mutation_extractor_manager import NUC_2_IDX, parse_line, quality_check


def test_split_reads():
    assert split_reads(".,A") == [".", ",", "A"]
    assert split_reads("^].$,+2AC.-1g*$") == ["^].$", ",+2AC", ".-1g", "*$"]
    assert split_reads("^+.^*,") == ["^+.", "^*,"]        # MAPQ chars that look like markers
    assert split_reads(".+12ACGTACGTACGTc") == [".+12ACGTACGTACGT", "c"]


def _cols(reads):
    """reads: (token, mapq, zn) -> annotated mpileup columns for one species."""
    return (str(len(reads)), "".join(t for t, _, _ in reads),
            "".join(chr(q + 33) for _, q, _ in reads), ",".join(str(z) for _, _, z in reads))


def test_zero_depth():
    assert site_intervals("0", "*", "*", "*") == (None, None)


def test_hand_worked():
    reads = [(".", 60, -1), (",", 20, -1), ("T", 5, -1), ("t", 3, 0)]
    plain, cont = site_intervals(*_cols(reads))
    # without continuity: '.' from t=6 (T at 5 conflicts) up to 60
    assert plain == (".", 6, 60)
    # with continuity: the rescued 't' is kept at every threshold and sets the base,
    # so the site is usable only once both '.' reads are gone (t > 60)
    assert cont == ("T", 61, ALWAYS_KEPT)


def test_indel_read_bounds_interval():
    reads = [(".", 60, -1), (".+1A", 30, -1)]
    assert site_intervals(*_cols(reads))[0] == (".", 31, 60)


def test_read_start_with_marker_mapq_is_not_an_indel():
    """'^' + MAPQ as a character is '*' at 9, '+' at 10 and '-' at 12. Since bf9d3c6 the read
    start is not read as an indel, so an agreeing one leaves the interval open."""
    for q in (9, 10, 12):
        reads = [(".", 60, -1), ("^" + chr(q + 33) + ".", q, -1)]
        assert site_intervals(*_cols(reads))[0] == (".", 0, 60)
        disagreeing = [(".", 60, -1), ("^" + chr(q + 33) + "T", q, -1)]
        assert site_intervals(*_cols(disagreeing))[0] == (".", q + 1, 60)


def _production_line(kept):
    """The pileup line a production run would see at this site: species 1 fixed and
    clean, species 2 holding only the kept reads (mpileup omits filtered reads)."""
    tokens = "".join(t for t, _, _ in kept)
    depth = len(kept)
    return f"chr1\t100\tA\t1\t.\tI\t{depth}\t{tokens or '*'}\t{'I' * depth or '*'}\n"


def _random_token(rng, mapq):
    tok = rng.choice(".,.,.,AaCcGgTt*")
    if rng.random() < 0.2:
        tok = "^" + chr(mapq + 33) + tok
    if rng.random() < 0.1:
        tok += rng.choice(["+1A", "-2cg", "+10ACGTACGTAC"])
    if rng.random() < 0.2:
        tok += "$"
    return tok


@pytest.mark.parametrize("seed", range(300))
def test_matches_quality_check_at_every_threshold(seed):
    rng = random.Random(seed)
    reads = []
    for _ in range(rng.randint(1, 5)):
        mapq = rng.choice([1, 3, 9, 10, 12, 20, 30, 60])
        reads.append((_random_token(rng, mapq), mapq, rng.choice([-1, -1, 0, 20, 60])))
    intervals = site_intervals(*_cols(reads))

    for continuity, interval in zip((False, True), intervals):
        for t in list(range(0, 62)) + [255]:
            kept = [r for r in reads if r[1] >= t or (continuity and r[2] >= 0)]
            line = parse_line(_production_line(kept))
            passes = quality_check(line, defaultdict(int))
            usable = interval is not None and interval[1] <= t <= interval[2]
            assert passes == usable, (reads, continuity, t, interval)
            if usable:
                base = line.nuc(NUC_2_IDX)
                assert interval[0] == ("." if base in ".," else base)


# ---------------------------------------------------------------- annotated scan

from coral.annotation import derive_counts, derive_run, scan_annotated_pileup, write_annotated_outputs
from coral.mutation_extractor_manager import scan_pileup
from coral.repeat_masker import RepeatMask


def _random_species_reads(rng, alt, carries_alt):
    """Mostly consistent reads for one species (its true base is alt or the
    reference), with occasional disagreeing bases, deletions and insertions."""
    true, other = (alt, ".") if carries_alt else (".", alt)
    reads = []
    for _ in range(rng.choice([0, 1, 1, 1, 2, 2, 3])):
        mapq = rng.choice([1, 9, 10, 20, 30, 60, 60])
        r = rng.random()
        tok = true if r < 0.9 else other if r < 0.96 else "*"
        if tok != "*" and rng.random() < 0.5:
            tok = "," if tok == "." else tok.lower()
        if rng.random() < 0.05:
            tok += "+1A"
        if rng.random() < 0.1:
            tok = "^" + chr(mapq + 33) + tok
        reads.append((tok, mapq, rng.choice([-1, -1, 0])))
    return reads


def _random_pileup(rng):
    """-> list of (chrom, pos, ref, reads_sp1, reads_sp2)."""
    sites = []
    for chrom in ("chr1", "chr2"):
        pos = 1
        for _ in range(rng.randint(20, 60)):
            pos += 1 if rng.random() < 0.85 else rng.randint(2, 3)
            ref = rng.choice("ACGT")
            alt = rng.choice([b for b in "ACGT" if b != ref])
            both = rng.random() < 0.05
            sites.append((chrom, pos, ref,
                          _random_species_reads(rng, alt, both or rng.random() < 0.05),
                          _random_species_reads(rng, alt, both or rng.random() < 0.05)))
    return sites


def _annotated_lines(sites):
    lines = []
    for chrom, pos, ref, r1, r2 in sites:
        if not r1 and not r2:
            continue                                   # mpileup omits empty positions
        cols = []
        for reads in (r1, r2):
            if reads:
                d, b, m, z = _cols(reads)
                cols += [d, b, "I" * len(reads), m, z]
            else:
                cols += ["0", "*", "*", "*", "*"]
        lines.append("\t".join([chrom, str(pos), ref] + cols) + "\n")
    return lines


def _production_lines(sites, t, continuity):
    lines = []
    for chrom, pos, ref, r1, r2 in sites:
        kept = [[r for r in reads if r[1] >= t or (continuity and r[2] >= 0)] for reads in (r1, r2)]
        if not kept[0] and not kept[1]:
            continue
        cols = []
        for k in kept:
            cols += [str(len(k)), "".join(tok for tok, _, _ in k) or "*", "I" * len(k) or "*"]
        lines.append("\t".join([chrom, str(pos), ref] + cols) + "\n")
    return lines


@pytest.mark.parametrize("seed", range(60))
@pytest.mark.parametrize("masked", [False, True])
def test_scan_matches_production_at_every_threshold(seed, masked):
    rng = random.Random(seed)
    sites = _random_pileup(rng)
    mask = RepeatMask({"chr1": [(4, 6)], "chr2": [(10, 10)]}) if masked else None

    calls = []
    counts = scan_annotated_pileup(iter(_annotated_lines(sites)),
                                   on_call=lambda *row: calls.append(row), ref_mask=mask)
    for mode in (0, 1):
        for t in (0, 1, 2, 9, 10, 11, 20, 21, 30, 31, 59, 60, 61, 255):
            # the mask only flags calls, so counts and calls equal an unmasked production scan
            rows1, rows2 = [], []
            mut1, mut2, trip1, trip2, summary, ref_diff = scan_pileup(
                iter(_production_lines(sites, t, mode)),
                on_mut1=lambda c, p, m: rows1.append((c, p, m)),
                on_mut2=lambda c, p, m: rows2.append((c, p, m)))
            expected = {"mut1": dict(mut1), "mut2": dict(mut2), "trip1": dict(trip1),
                        "trip2": dict(trip2), "classes": summary["site_classes"],
                        "ref_diff": dict(ref_diff)}
            assert derive_counts(counts, t, mode) == expected, (seed, mode, t)

            def derived_rows(keep_repeats):
                return {taxon: [(c, p, m) for c, p, md, tx, m, lo, hi, rep in calls
                                if md == mode and tx == taxon and lo <= t <= hi and (keep_repeats or not rep)]
                        for taxon in (1, 2)}
            rows = derived_rows(keep_repeats=True)
            assert rows[1] == rows1 and rows[2] == rows2, (seed, mode, t)

            if mask:
                # calls not flagged as repeats are exactly production's masked calls
                masked1, masked2 = [], []
                scan_pileup(iter(_production_lines(sites, t, mode)),
                            on_mut1=lambda c, p, m: masked1.append((c, p, m)),
                            on_mut2=lambda c, p, m: masked2.append((c, p, m)), ref_mask=mask)
                rows = derived_rows(keep_repeats=False)
                assert rows[1] == masked1 and rows[2] == masked2, (seed, mode, t)


def test_repeat_flag_marks_calls_whose_window_touches_the_mask():
    flagged = unflagged = 0
    for seed in range(60):
        sites = _random_pileup(random.Random(seed))
        mask = RepeatMask({"chr1": [(4, 6)], "chr2": [(10, 10)]})
        calls = []
        scan_annotated_pileup(iter(_annotated_lines(sites)), on_call=lambda *r: calls.append(r), ref_mask=mask)
        for c, p, *_, rep in calls:
            assert rep == int(any(mask.contains(c, q) for q in (p - 1, p, p + 1)))
            flagged += rep
            unflagged += 1 - rep
        no_mask = []
        scan_annotated_pileup(iter(_annotated_lines(sites)), on_call=lambda *r: no_mask.append(r))
        assert all(r[-1] == 0 for r in no_mask)
    assert flagged > 0 and unflagged > 0


def test_scan_random_cases_are_not_trivial():
    """The randomized comparison is only meaningful if production actually calls
    mutations and classifies sites in those cases."""
    classes = defaultdict(int)
    for seed in range(60):
        sites = _random_pileup(random.Random(seed))
        for (mode, site_class, label, lo, hi), n in scan_annotated_pileup(iter(_annotated_lines(sites))).items():
            classes[site_class] += n
    assert classes["taxa1_mut"] > 5 and classes["taxa2_mut"] > 5 and classes["identical"] > 50
    assert classes["ref_differs"] > 0 and classes["flanks_not_conserved"] > 0


def test_write_annotated_outputs_round_trip(tmp_path):
    import gzip
    from coral.annotation import CALLS_HEADER, COUNTS_HEADER, write_annotated_outputs

    lines = _annotated_lines(_random_pileup(random.Random(17)))   # seed with several calls
    pileup = tmp_path / "x.annotated.pileup.gz"
    with gzip.open(pileup, "wt") as f:
        f.writelines(lines)

    calls = []
    counts = scan_annotated_pileup(iter(lines), on_call=lambda *r: calls.append(r))
    assert calls, "test data should contain calls"

    calls_path, counts_path = write_annotated_outputs(str(pileup), str(tmp_path / "Annotated"),
                                                      "ref", "sp1", "sp2", verbose=False)
    with gzip.open(calls_path, "rt") as f:
        assert f.readline() == CALLS_HEADER
        rows = [line.rstrip("\n").split(",") for line in f]
    names = {1: "sp1", 2: "sp2"}
    assert rows == [[c, str(p), str(m), names[t], mut, str(lo), str(hi), str(rep)]
                    for c, p, m, t, mut, lo, hi, rep in calls]

    with gzip.open(counts_path, "rt") as f:
        assert f.readline() == COUNTS_HEADER
        read_back = {}
        for line in f:
            mode, site_class, label, lo, hi, n = line.rstrip("\n").split("\t")
            read_back[(int(mode), site_class, label, int(lo), int(hi))] = int(n)
    assert read_back == dict(counts)

    # cached: a second call must not rewrite
    mtime = (tmp_path / "Annotated" / "ref__sp1__sp2__calls.csv.gz").stat().st_mtime_ns
    write_annotated_outputs(str(pileup), str(tmp_path / "Annotated"), "ref", "sp1", "sp2", verbose=False)
    assert (tmp_path / "Annotated" / "ref__sp1__sp2__calls.csv.gz").stat().st_mtime_ns == mtime


def _synthetic_run(tmp_path, seed, low_mapq=1):
    """A minimal finished --annotate run directory built from random sites."""
    import gzip, json
    rng = random.Random(seed)
    sites = [(f"s{i}_{chrom}", pos, ref, r1, r2)                  # 30 pileups, separate contigs
             for i in range(30) for chrom, pos, ref, r1, r2 in _random_pileup(rng)]
    run_dir = tmp_path / "run"
    (run_dir).mkdir()
    pileup = run_dir / "x.annotated.pileup.gz"
    with gzip.open(pileup, "wt") as f:
        f.writelines(_annotated_lines(sites))
    write_annotated_outputs(str(pileup), str(run_dir / "Annotated"), "ref", "sp1", "sp2", verbose=False)
    with open(run_dir / "run_summary.json", "w") as f:
        json.dump({"run": {"outgroup": "ref", "taxa": ["sp1", "sp2"],
                           "parameters": {"annotate": True, "low_mapq": low_mapq}}}, f)
    return sites, run_dir


@pytest.mark.parametrize("seed", [17, 3, 42])
@pytest.mark.parametrize("mapq,continuity", [(1, False), (20, False), (30, True), (60, True), (60, False)])
def test_derive_run_writes_what_production_writes(tmp_path, seed, mapq, continuity):
    import gzip, json
    import pandas as pd
    from coral.mutation_extractor_manager import MutationExtractor, MutationNormalizer

    sites, run_dir = _synthetic_run(tmp_path, seed)
    derived = derive_run(str(run_dir), mapq, continuity, verbose=False)

    prod = tmp_path / "production"
    prod.mkdir()
    with gzip.open(prod / "p.pileup.gz", "wt") as f:
        f.writelines(_production_lines(sites, mapq, continuity))
    MutationExtractor("ref", "sp1", "sp2", str(prod / "p.pileup.gz"), str(prod / "Mutations"),
                      str(prod / "Triplets"), verbose=False).extract()
    MutationNormalizer(input_dir=str(prod), output_dir=str(prod / "Tables"), verbose=False).normalize()

    for sub, names in [("Mutations", ["sp1__sp2__ref__mutations.json", "sp2__sp1__ref__mutations.json"]),
                       ("Triplets", ["sp1__sp2__ref__triplets.json", "sp2__sp1__ref__triplets.json"])]:
        for name in names:
            assert json.load(open(f"{derived}/{sub}/{name}")) == json.load(open(prod / sub / name)), name
    for name in ["sp1__sp2__ref__mutations.csv.gz", "sp2__sp1__ref__mutations.csv.gz"]:
        assert gzip.open(f"{derived}/Mutations/{name}").read() == gzip.open(prod / "Mutations" / name).read()

    d_summary, p_summary = json.load(open(f"{derived}/run_summary.json")), json.load(open(prod / "run_summary.json"))
    assert d_summary["site_classes"] == p_summary["site_classes"]
    assert d_summary["reference_difference_spectrum"] == p_summary["reference_difference_spectrum"]
    assert d_summary["derived_from"]["mapq"] == mapq

    for table in ["collapsed_mutations.tsv", "normalized_scaled.tsv", "scaled_raw.tsv", "triplets.tsv"]:
        a = pd.read_csv(f"{derived}/Tables/{table}", sep="\t", index_col=0)
        b = pd.read_csv(prod / "Tables" / table, sep="\t", index_col=0)
        pd.testing.assert_frame_equal(a.sort_index().sort_index(axis=1), b.sort_index().sort_index(axis=1))


def test_derive_run_refuses_below_low_mapq_and_non_annotate_runs(tmp_path):
    import json
    _, run_dir = _synthetic_run(tmp_path, 17, low_mapq=5)
    with pytest.raises(ValueError, match="below"):
        derive_run(str(run_dir), 4, True, verbose=False)
    summary = json.load(open(run_dir / "run_summary.json"))
    summary["run"]["parameters"]["annotate"] = False
    json.dump(summary, open(run_dir / "run_summary.json", "w"))
    with pytest.raises(ValueError, match="not an --annotate run"):
        derive_run(str(run_dir), 30, True, verbose=False)


# ---------------------------------------------------------------- pair mode

from coral.annotation import _change_triplets, derive_pair_counts
from coral.mutation_extractor_manager import PairExtractor, fold, scan_pair


def _pair_sites(rng):
    """One target species against the reference: (chrom, pos, ref, reads)."""
    return [(chrom, pos, ref, reads) for chrom, pos, ref, reads, _ in _random_pileup(rng)]


def _annotated_pair_lines(sites):
    """A one-sample annotated pileup: mpileup writes no line where nothing aligns."""
    return ["\t".join([chrom, str(pos), ref] + [d, b, "I" * len(reads), m, z]) + "\n"
            for chrom, pos, ref, reads in sites if reads
            for d, b, m, z in [_cols(reads)]]


def _production_pair_lines(sites, t, continuity):
    lines = []
    for chrom, pos, ref, reads in sites:
        kept = [r for r in reads if r[1] >= t or (continuity and r[2] >= 0)]
        if not kept:
            continue
        lines.append("\t".join([chrom, str(pos), ref, str(len(kept)),
                                "".join(tok for tok, _, _ in kept), "I" * len(kept)]) + "\n")
    return lines


@pytest.mark.parametrize("seed", range(60))
def test_pair_scan_matches_production_at_every_threshold(seed):
    sites = _pair_sites(random.Random(seed))
    calls = []
    counts = scan_annotated_pileup(iter(_annotated_pair_lines(sites)), n_species=1,
                                   on_call=lambda *row: calls.append(row))
    for mode in (0, 1):
        for t in (0, 1, 2, 9, 10, 11, 20, 21, 30, 31, 59, 60, 61, 255):
            rows = []
            muts, trips, summary = scan_pair(iter(_production_pair_lines(sites, t, mode)),
                                             on_mut=lambda c, p, change, m: rows.append((c, p, change, m)))
            assert derive_pair_counts(counts, t, mode) == {
                "mut": dict(muts), "trip": dict(trips), "classes": summary["site_classes"]}, (seed, mode, t)
            derived = [(c, p, change, fold(*_change_triplets(change)))
                       for c, p, md, _, change, lo, hi, _ in calls if md == mode and lo <= t <= hi]
            assert derived == rows, (seed, mode, t)


def test_pair_scan_random_cases_are_not_trivial():
    classes = defaultdict(int)
    for seed in range(60):
        sites = _pair_sites(random.Random(seed))
        for (mode, site_class, label, lo, hi), n in scan_annotated_pileup(
                iter(_annotated_pair_lines(sites)), n_species=1).items():
            classes[site_class] += n
    assert classes["differs"] > 5 and classes["identical"] > 50 and classes["flanks_not_conserved"] > 0


def test_pair_scan_classifies_non_acgt_like_production():
    """Agreeing N reads are clean for sample_problem, and both scans call the window not_acgt."""
    sites = [("chr1", p, "ACGTA"[p - 1], [("N" if p == 3 else ".", 60, -1)]) for p in range(1, 6)]
    counts = scan_annotated_pileup(iter(_annotated_pair_lines(sites)), n_species=1)
    muts, trips, summary = scan_pair(iter(_production_pair_lines(sites, 60, False)))
    assert summary["site_classes"]["not_acgt"] == 3
    assert derive_pair_counts(counts, 60, False)["classes"] == summary["site_classes"]


def _synthetic_pair_run(tmp_path, seed, low_mapq=1):
    """A minimal finished pair --annotate run directory."""
    import gzip, json
    rng = random.Random(seed)
    sites = [(f"s{i}_{chrom}", pos, ref, reads)
             for i in range(30) for chrom, pos, ref, reads in _pair_sites(rng)]
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    pileup = run_dir / "x.annotated.pileup.gz"
    with gzip.open(pileup, "wt") as f:
        f.writelines(_annotated_pair_lines(sites))
    write_annotated_outputs(str(pileup), str(run_dir / "Annotated"), "ref", "sp1", verbose=False)
    with open(run_dir / "run_summary.json", "w") as f:
        json.dump({"run": {"outgroup": "ref", "taxa": ["sp1"],
                           "parameters": {"annotate": True, "low_mapq": low_mapq}}}, f)
    return sites, run_dir


@pytest.mark.parametrize("seed", [17, 3, 42])
@pytest.mark.parametrize("mapq,continuity", [(1, False), (20, False), (30, True), (60, True), (60, False)])
def test_pair_derive_run_writes_what_production_writes(tmp_path, seed, mapq, continuity):
    import gzip, json

    sites, run_dir = _synthetic_pair_run(tmp_path, seed)
    derived = derive_run(str(run_dir), mapq, continuity, verbose=False)

    prod = tmp_path / "production"
    prod.mkdir()
    with gzip.open(prod / "p.pileup.gz", "wt") as f:
        f.writelines(_production_pair_lines(sites, mapq, continuity))
    PairExtractor("ref", "sp1", None, None, str(prod / "p.pileup.gz"),
                  str(prod / "Mutations"), str(prod / "Triplets"), verbose=False).extract()

    for sub, name in [("Mutations", "sp1__ref__mutations.json"), ("Triplets", "sp1__ref__triplets.json")]:
        assert json.load(open(f"{derived}/{sub}/{name}")) == json.load(open(prod / sub / name)), name
    assert (gzip.open(f"{derived}/Mutations/sp1__ref__mutations.csv.gz").read()
            == gzip.open(prod / "Mutations" / "sp1__ref__mutations.csv.gz").read())

    d_summary = json.load(open(f"{derived}/run_summary.json"))
    p_summary = json.load(open(prod / "run_summary.json"))
    assert d_summary["site_classes"] == p_summary["site_classes"]
    assert d_summary["derived_from"]["mapq"] == mapq
    assert not os.path.exists(f"{derived}/Tables")      # a pair run writes no Tables/


def test_pair_annotate_takes_the_pileup_path(tmp_path):
    """--annotate on a pair scans the annotated pileup; the BAM scan (extract_pair) is skipped."""
    from coral.pipeline import MutationExtractionPipeline

    pipeline = MutationExtractionPipeline(species_list=[("sp1", "ACC1")], outgroup=("ref", "ACC0"),
                                          base_output_dir=str(tmp_path), annotate=True, verbose=False)
    assert pipeline.pair
    stages = []
    for name in ("download_genomes", "download_index_and_fragment_genomes", "align_species", "generate_pileup",
                 "extract_annotated", "extract_pair", "extract_mutations_and_triplets",
                 "extract_intervals", "run_plots", "cleanup", "_collect_genome_stats",
                 "_write_run_summary"):
        setattr(pipeline, name, (lambda n: lambda *a, **k: stages.append(n))(name))
    os.makedirs(pipeline.output_dir, exist_ok=True)
    pipeline.run()

    assert "generate_pileup" in stages and "extract_annotated" in stages
    assert "extract_pair" not in stages and "extract_mutations_and_triplets" not in stages


@pytest.mark.parametrize("seed", range(20))
@pytest.mark.parametrize("k", [2, 5])
def test_pair_indel_window_matches_production(seed, k):
    """With --indel-window k the annotated scan flags each window with the highest MAPQ of a
    nearby indel read; dropping the windows whose flag reaches the threshold must leave exactly
    what a production scan at that threshold keeps."""
    sites = _pair_sites(random.Random(seed))
    counts = scan_annotated_pileup(iter(_annotated_pair_lines(sites)), n_species=1, indel_window=k)
    for mode in (0, 1):
        for t in (1, 10, 30, 60):
            muts, trips, summary = scan_pair(iter(_production_pair_lines(sites, t, mode)), indel_window=k)
            assert derive_pair_counts(counts, t, mode, indel=True) == {
                "mut": dict(muts), "trip": dict(trips), "classes": summary["site_classes"]}, (seed, k, mode, t)


def test_no_cleanup_keeps_the_intermediates(tmp_path):
    """--no-cleanup skips the stage that deletes the raw BAMs, the pileup and the genome folders."""
    from coral.pipeline import MutationExtractionPipeline

    def stages_of(**params):
        pipeline = MutationExtractionPipeline(species_list=[("sp1", "ACC1")], outgroup=("ref", "ACC0"),
                                              base_output_dir=str(tmp_path), verbose=False, **params)
        seen = []
        for name in ("download_genomes", "download_index_and_fragment_genomes", "align_species", "generate_pileup",
                     "extract_annotated", "extract_pair", "extract_intervals", "cleanup",
                     "_collect_genome_stats", "_write_run_summary"):
            setattr(pipeline, name, (lambda n: lambda *a, **k: seen.append(n))(name))
        os.makedirs(pipeline.output_dir, exist_ok=True)
        pipeline.run()
        return seen

    assert "cleanup" in stages_of()
    assert "cleanup" not in stages_of(no_cleanup=True)


def test_parallel_annotate_matches_the_whole_genome_scan(tmp_path, monkeypatch):
    """A chromosome per process must give the file a single pass over the whole-genome
    pileup gives, byte for byte -- including the order the chromosomes are written in."""
    import gzip
    from coral import annotation

    chroms = ["chrA", "chrB", "chrC"]
    lines, paths = {}, {}
    for i, c in enumerate(chroms):
        sites = [(c, pos, ref, reads)
                 for ch, pos, ref, reads in _pair_sites(random.Random(i)) if ch == "chr1"]
        lines[c] = _annotated_pair_lines(sites)
        paths[c] = str(tmp_path / f"{c}.pileup")
        with open(paths[c], "w") as f:
            f.writelines(lines[c])

    whole = tmp_path / "whole.annotated.pileup.gz"
    with gzip.open(whole, "wt") as f:                     # chromosomes in .fai order
        for c in chroms:
            f.writelines(lines[c])
    fai = tmp_path / "ref.fasta.fai"
    fai.write_text("".join(f"{c}\t{1000 * (3 - i)}\t0\t60\t61\n" for i, c in enumerate(chroms)))

    monkeypatch.setattr(annotation, "mpileup_cmd",
                        lambda ref_fasta, bams, region=None, annotate=False: ["cat", paths[region]])
    monkeypatch.setattr(annotation, "_chroms_with_reads", lambda bams: set(chroms))

    serial = write_annotated_outputs(str(whole), str(tmp_path / "serial"), "ref", "sp1", verbose=False)
    parallel = write_annotated_outputs(str(whole), str(tmp_path / "parallel"), "ref", "sp1", verbose=False,
                                       ref_fasta=str(tmp_path / "ref.fasta"), bams=["fake.bam"], cores=3)
    for a, b in zip(serial, parallel):
        assert gzip.open(a).read() == gzip.open(b).read(), os.path.basename(a)
    assert b"chrA" in gzip.open(parallel[0]).read()       # the scan found calls to write
