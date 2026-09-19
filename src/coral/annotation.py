"""Annotate mode: record, for every site, the MAPQ thresholds at which it is callable.

A read is kept at threshold t if MAPQ >= t, or with continuity if ZN >= 0.
A site is callable when its kept reads agree on one base and none has '*' or '+'.
With an indel window k > 1, each window also records indel_mapq: the highest MAPQ (256 when rescued by
continuity) of a read with '*' or '+' within k bp of the middle base.
"""
import gzip
import json
import os
from collections import defaultdict

from .mutation_extractor_manager import MutationNormalizer, _write_extractor_outputs, detect_mutation_triplet, scan_windows
from .utils import log

ALWAYS_KEPT = 256       # MAPQ given to continuity-rescued reads
CALLS_HEADER = "chromosome,position,continuity,taxon,mutation,mapq_lo,mapq_hi,repeat\n"
COUNTS_HEADER = "continuity\tsite_class\tlabel\tmapq_lo\tmapq_hi\tcount\n"
CALLS_HEADER_INDEL = CALLS_HEADER[:-1] + ",indel_mapq\n"                      # with an indel window k > 1
COUNTS_HEADER_INDEL = "continuity\tsite_class\tlabel\tmapq_lo\tmapq_hi\tindel_mapq\tcount\n"
SPECIES_COLS = ((3, 4, 6, 7), (8, 9, 11, 12))   # depth, bases, MAPQ, ZN columns per species


def split_reads(bases):
    """Split an mpileup bases field into one token per read."""
    tokens, i = [], 0
    while i < len(bases):
        start = i
        i += 3 if bases[i] == '^' else 1
        if i < len(bases) and bases[i] in '+-':
            j = i + 1
            while bases[j].isdigit():
                j += 1
            i = j + int(bases[i + 1:j])
        if i < len(bases) and bases[i] == '$':
            i += 1
        tokens.append(bases[start:i])
    return tokens


def _interval(reads):
    """(mapq, base, rejects) reads -> (base, lo, hi) callable threshold range, or None."""
    if not reads:
        return None
    hi, base, _ = max(reads)
    lo = max([q for q, b, rejects in reads if rejects or b != base], default=-1) + 1
    return (base, lo, hi) if lo <= hi else None


def site_intervals(depth, bases, mapqs, zns):
    """-> (interval without continuity, interval with continuity), each (base, lo, hi) or None."""
    if depth == "0":
        return None, None
    plain, rescued = [], []
    for token, q, zn in zip(split_reads(bases), mapqs, zns.split(',')):
        base = token[2] if token[0] == '^' else token[0]
        base = '.' if base in '.,' else base.upper()
        body = token[2:] if token[0] == '^' else token   # '^' + MAPQ as a character, which can be * or +
        rejects = '*' in body or '+' in body
        plain.append((ord(q) - 33, base, rejects))
        rescued.append((ALWAYS_KEPT if int(zn) >= 0 else ord(q) - 33, base, rejects))
    return _interval(plain), _interval(rescued)


def _parse_annotated(line):
    parts = line.rstrip('\n').split('\t')
    if len(parts) != 13:
        return None
    return parts[0], int(parts[1]), parts[2][0].upper(), [site_intervals(*(parts[c] for c in cols)) for cols in SPECIES_COLS]


def _indel_mapqs(line, parsed=None):
    """-> (highest MAPQ, highest MAPQ with continuity rescue) of the reads with '*' or '+' at this base, or None."""
    parts = line.rstrip('\n').split('\t')
    plain = rescued = -1
    for depth, bases, mapqs, zns in ((parts[c] for c in cols) for cols in SPECIES_COLS):
        if depth == '0' or ('*' not in bases and '+' not in bases):   # no reads: a lone '*', not a deletion
            continue
        for token, q, zn in zip(split_reads(bases), mapqs, zns.split(',')):
            body = token[2:] if token[0] == '^' else token
            if '*' in body or '+' in body:
                plain = max(plain, ord(q) - 33)
                rescued = max(rescued, ALWAYS_KEPT if int(zn) >= 0 else ord(q) - 33)
    return (plain, rescued) if plain >= 0 else None


def _window(lines, mode):
    """Intersect the 3 positions x 2 species intervals -> (triplets, lo, hi) or None."""
    lo, hi, triplets = 0, ALWAYS_KEPT, [[], [], []]
    for _, _, ref, species in lines:
        triplets[0].append(ref)
        for k, intervals in enumerate(species, start=1):
            if intervals[mode] is None:
                return None
            base, s_lo, s_hi = intervals[mode]
            lo, hi = max(lo, s_lo), min(hi, s_hi)
            triplets[k].append(ref if base == '.' else base)
    return (triplets, lo, hi) if lo <= hi else None


def scan_annotated_pileup(line_iter, on_call=None, ref_mask=None, indel_window=1):
    """Classify every 3-position window for continuity off (0) and on (1).

    Returns counts keyed (mode, site_class, label, lo, hi), plus indel_mapq when indel_window > 1.
    Calls go to on_call(chrom, pos, mode, taxon, mutation, lo, hi, repeat[, indel_mapq]), where
    repeat is 1 if any window position is in ref_mask.
    """
    counts = defaultdict(int)
    both_max = lambda a, b: (max(a[0], b[0]), max(a[1], b[1]))
    for window, near in scan_windows(line_iter, _parse_annotated, lambda line: line is not None,
                                     _indel_mapqs, indel_window, combine=both_max):
        (c0, p0, _, _), (c1, p1, _, _), (c2, p2, _, _) = window
        for mode in (0, 1):
            found = _window(window, mode)
            if found is None:
                continue
            triplets, lo, hi = found
            m1, m2, _, _, site_class = detect_mutation_triplet(triplets)
            ref = triplets[0]
            label = {'taxa1_mut': m1, 'taxa2_mut': m2, 'identical': ''.join(ref),
                     'ref_differs': f"{ref[0]}[{triplets[1][1]}>{ref[1]}]{ref[2]}"}.get(site_class, '')
            if indel_window == 1:
                flag = ()
                counts[(mode, site_class, label, lo, hi)] += 1
            else:
                flag = (near[mode] if near else 0,)
                counts[(mode, site_class, label, lo, hi, flag[0])] += 1
            if on_call and (m1 or m2):
                repeat = int(bool(ref_mask) and any(ref_mask.contains(c1, p) for p in (p0, p1, p2)))
                on_call(c1, p1, mode, 1 if m1 else 2, m1 or m2, lo, hi, repeat, *flag)
    return counts


def write_annotated_outputs(pileup_path, output_dir, reference, taxon1, taxon2,
                            ref_mask=None, no_cache=False, verbose=True, indel_window=1):
    """Scan an annotated pileup into <reference>__<taxon1>__<taxon2>__calls.csv.gz
    (one row per call and continuity mode) and __site_counts.tsv.gz."""
    os.makedirs(output_dir, exist_ok=True)
    stem = os.path.join(output_dir, f"{reference}__{taxon1}__{taxon2}")
    calls_path, counts_path = f"{stem}__calls.csv.gz", f"{stem}__site_counts.tsv.gz"
    if not no_cache and os.path.exists(calls_path) and os.path.exists(counts_path):
        log(f"Annotated outputs already exist: {calls_path}", verbose)
        return calls_path, counts_path

    taxa = (None, taxon1, taxon2)
    with gzip.open(pileup_path, 'rt') as pileup, gzip.open(calls_path + ".tmp", 'wt') as calls:
        calls.write(CALLS_HEADER if indel_window == 1 else CALLS_HEADER_INDEL)
        counts = scan_annotated_pileup(
            pileup, ref_mask=ref_mask, indel_window=indel_window,
            on_call=lambda chrom, pos, mode, taxon, mutation, lo, hi, repeat, *flag:
                calls.write(",".join(map(str, (chrom, pos, mode, taxa[taxon], mutation, lo, hi, repeat) + flag)) + "\n"))
    with gzip.open(counts_path + ".tmp", 'wt') as out:
        out.write(COUNTS_HEADER if indel_window == 1 else COUNTS_HEADER_INDEL)
        for key in sorted(counts):
            out.write("\t".join(map(str, key + (counts[key],))) + "\n")
    os.replace(calls_path + ".tmp", calls_path)
    os.replace(counts_path + ".tmp", counts_path)
    log(f"Annotated outputs written to: {output_dir}", verbose)
    return calls_path, counts_path


def read_site_counts(path):
    """-> {(mode, site_class, label, lo, hi[, indel_mapq]): count}"""
    with gzip.open(path, 'rt') as f:
        header = f.readline()
        if header not in (COUNTS_HEADER, COUNTS_HEADER_INDEL):
            raise ValueError(f"unexpected header in {path}")
        counts = {}
        for line in f:
            mode, site_class, label, lo, hi, *flag, n = line.rstrip('\n').split('\t')
            counts[(int(mode), site_class, label, int(lo), int(hi), *map(int, flag))] = int(n)
    return counts


def read_calls(path):
    """-> iterator of (chrom, pos, mode, taxon, mutation, lo, hi, repeat), plus indel_mapq when the file has it"""
    with gzip.open(path, 'rt') as f:
        header = f.readline()
        if header not in (CALLS_HEADER, CALLS_HEADER_INDEL):
            raise ValueError(f"unexpected header in {path}")
        for line in f:
            chrom, pos, mode, taxon, mutation, lo, hi, repeat, *flag = line.rstrip('\n').split(',')
            yield (chrom, int(pos), int(mode), taxon, mutation, int(lo), int(hi), int(repeat), *map(int, flag))


def derive_counts(counts, mapq, continuity, indel=False):
    """Scan outputs at one setting: mut1, mut2, trip1, trip2, classes and ref_diff.
    With indel=True, windows flagged as near an indel at this threshold are left out."""
    mode = int(bool(continuity))
    out = {name: defaultdict(int) for name in ("mut1", "mut2", "trip1", "trip2", "classes", "ref_diff")}
    for (m, site_class, label, lo, hi, *flag), n in counts.items():
        if m != mode or not lo <= mapq <= hi or (indel and flag and flag[0] >= mapq):
            continue
        out["classes"][site_class] += n
        if site_class in ("taxa1_mut", "taxa2_mut", "identical"):
            if site_class != "identical":
                out["mut1" if site_class == "taxa1_mut" else "mut2"][label] += n
            context = label if site_class == "identical" else label[0] + label[2] + label[-1]
            out["trip1"][context] += n
            out["trip2"][context] += n
        elif site_class == "ref_differs":
            out["ref_diff"][label] += n
    return {name: dict(d) for name, d in out.items()}


def derive_run(run_dir, mapq, continuity, output_dir=None, no_cache=False, verbose=True, indel=False):
    """Write the Mutations/, Triplets/, run_summary.json and Tables/ that a normal run at
    this setting would produce, into run_dir/Derived/mapq<T>_<continuity|no_continuity>[_indel].
    mapq must be at least the run's --low-mapq. indel=True applies the run's indel window flag."""
    with open(os.path.join(run_dir, "run_summary.json")) as f:
        run = json.load(f)["run"]
    params = run.get("parameters", {})
    if not params.get("annotate"):
        raise ValueError(f"{run_dir} is not an --annotate run")
    if mapq < params.get("low_mapq", 1):
        raise ValueError(f"mapq {mapq} is below this run's low_mapq {params.get('low_mapq', 1)}")
    if indel and params.get("indel_window", 1) <= 1:
        raise ValueError(f"{run_dir} was annotated without an indel window (--indel-window > 1)")

    reference, (taxon1, taxon2) = run["outgroup"], run["taxa"]
    setting = f"mapq{mapq}_{'continuity' if continuity else 'no_continuity'}" + ("_indel" if indel else "")
    output_dir = output_dir or os.path.join(run_dir, "Derived", setting)
    summary_json = os.path.join(output_dir, "run_summary.json")
    if os.path.exists(summary_json) and not no_cache:
        log(f"Derived outputs already exist: {output_dir}", verbose)
        return output_dir

    stem = os.path.join(run_dir, "Annotated", f"{reference}__{taxon1}__{taxon2}")
    derived = derive_counts(read_site_counts(f"{stem}__site_counts.tsv.gz"), mapq, continuity, indel)
    mut_dir, trip_dir = os.path.join(output_dir, "Mutations"), os.path.join(output_dir, "Triplets")
    os.makedirs(mut_dir, exist_ok=True)
    os.makedirs(trip_dir, exist_ok=True)
    pair = {taxon1: f"{taxon1}__{taxon2}__{reference}", taxon2: f"{taxon2}__{taxon1}__{reference}"}

    mode = int(bool(continuity))
    csvs = {taxon: gzip.open(os.path.join(mut_dir, f"{name}__mutations.csv.gz"), 'wt') for taxon, name in pair.items()}
    try:
        for csv in csvs.values():
            csv.write("chromosome,position,mutation\n")
        for chrom, pos, m, taxon, mutation, lo, hi, _, *near in read_calls(f"{stem}__calls.csv.gz"):
            if m == mode and lo <= mapq <= hi and not (indel and near and near[0] >= mapq):
                csvs[taxon].write(f"{chrom},{pos},{mutation}\n")
    finally:
        for csv in csvs.values():
            csv.close()

    # run_summary.json is written last: it marks a complete derivation
    _write_extractor_outputs(
        derived["mut1"], derived["mut2"], derived["trip1"], derived["trip2"],
        os.path.join(mut_dir, f"{pair[taxon1]}__mutations.json"),
        os.path.join(mut_dir, f"{pair[taxon2]}__mutations.json"),
        os.path.join(trip_dir, f"{pair[taxon1]}__triplets.json"),
        os.path.join(trip_dir, f"{pair[taxon2]}__triplets.json"),
        summary={"site_classes": derived["classes"],
                 "derived_from": {"run_dir": os.path.abspath(run_dir), "mapq": mapq, "continuity": bool(continuity),
                                  **({"indel": True} if indel else {})}},
        summary_json=summary_json + ".tmp", ref_diff=derived["ref_diff"])
    MutationNormalizer(input_dir=output_dir, output_dir=os.path.join(output_dir, "Tables"),
                       divergence_time=params.get("divergence_time"), verbose=verbose).normalize()
    os.replace(summary_json + ".tmp", summary_json)
    log(f"Derived {setting} outputs written to: {output_dir}", verbose)
    return output_dir
