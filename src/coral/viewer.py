"""Self-contained HTML browser for a finished `coral run_single` run.

Along the outgroup (reference) coordinates the page shows the reference repeat mask, the
reference sequence and each ingroup branch's mutation calls. A call is flagged when every
read behind it came from sequence masked in that branch's own genome, and given a second
run made without masking (`compare`), the calls the mask removed are shown as well.

It takes two steps because the pipeline's cleanup deletes the genomes:
  * `preserve_viewer_inputs` (run_single --viewer, just before cleanup) keeps the
    reference FASTA (bgzipped, faidx-indexed) and every repeat-mask BED under
    <run>/Viewer/inputs, with a manifest.
  * `build_viewer` (`coral view`, or automatically at the end of a --viewer run) reads
    those, the final BAMs and Mutations/, and writes one HTML file with the data embedded.
"""
import base64
import csv
import gzip
import json
import os
import re
import shutil
from importlib import resources

import numpy as np
import pysam

from .repeat_masker import RepeatMask
from .utils import log

VIEWER_DIR = "Viewer"
INPUTS_DIR = os.path.join(VIEWER_DIR, "inputs")
MANIFEST = "manifest.json"
TEMPLATE = "viewer_template.html"
DEFAULT_MAX_GENOME_BP = 50_000_000     # ~17 MB of embedded sequence; larger genomes use --region


class ViewerError(Exception):
    """The run can't be turned into a page as asked: inputs missing, too large, bad region."""


# ---- step 1: keep what cleanup would delete ----

def preserve_viewer_inputs(output_dir, reference, genomes, mask_beds, parameters=None, verbose=True):
    """Copy the reference FASTA (bgzipped + faidx) and each genome's repeat-mask BED into
    <output_dir>/Viewer/inputs and write a manifest. `reference` and `genomes` need
    `.name`, `.accession` and (reference only) `.fasta_path`; `mask_beds` maps genome
    name -> BED path, or None for genomes that were not masked."""
    inputs = os.path.join(output_dir, INPUTS_DIR)
    os.makedirs(inputs, exist_ok=True)

    fasta = os.path.join(inputs, "reference.fa.gz")
    pysam.tabix_compress(reference.fasta_path, fasta, force=True)
    for ext in (".fai", ".gzi"):
        if os.path.exists(fasta + ext):
            os.remove(fasta + ext)
    pysam.faidx(fasta)

    masks = {}
    for name, bed in (mask_beds or {}).items():
        if bed and os.path.exists(bed):
            dest = f"mask.{name}.bed"
            shutil.copyfile(bed, os.path.join(inputs, dest))
            masks[name] = dest

    manifest = {
        "reference": reference.name,
        "taxa": [g.name for g in genomes],
        "accessions": {g.name: g.accession for g in [reference, *genomes]},
        "fasta": os.path.basename(fasta),
        "masks": masks,
        "parameters": {k: v for k, v in (parameters or {}).items()
                       if v is None or isinstance(v, (str, int, float, bool))},
    }
    with open(os.path.join(inputs, MANIFEST), "w") as fh:
        json.dump(manifest, fh, indent=2)
    log(f"Viewer inputs saved to: {inputs}", verbose)
    return inputs


# ---- encoding helpers ----

_LUT = np.full(256, 4, np.uint8)
for _i, _b in enumerate("ACGT"):
    _LUT[ord(_b)] = _LUT[ord(_b.lower())] = _i


def pack_sequence(seq):
    """2-bit pack `seq` (A0 C1 G2 T3, most-significant pair first) as base64. Non-ACGT
    bases are returned separately as flat [offset, length, ...] runs (0-based offsets)."""
    codes = _LUT[np.frombuffer(seq.encode("ascii", "replace"), np.uint8)]
    other = codes == 4
    runs = []
    if other.any():
        edges = np.diff(np.r_[0, other.view(np.int8), 0])
        starts, ends = np.flatnonzero(edges == 1), np.flatnonzero(edges == -1)
        runs = np.column_stack([starts, ends - starts]).ravel().tolist()
    codes = np.where(other, 0, codes).astype(np.uint8)
    codes = np.pad(codes, (0, -len(codes) % 4)).reshape(-1, 4)
    packed = (codes[:, 0] << 6) | (codes[:, 1] << 4) | (codes[:, 2] << 2) | codes[:, 3]
    return base64.b64encode(packed.astype(np.uint8).tobytes()).decode("ascii"), runs


def unpack_sequence(b64, length, runs=()):
    """Inverse of `pack_sequence` (runs as 0-based offsets), for tests and debugging."""
    packed = np.frombuffer(base64.b64decode(b64), np.uint8)
    codes = np.stack([(packed >> s) & 3 for s in (6, 4, 2, 0)], axis=1).ravel()[:length]
    seq = np.array(list("ACGT"))[codes]
    for off, n in zip(runs[0::2], runs[1::2]):
        seq[off:off + n] = "N"
    return "".join(seq)


_REGION = re.compile(r"^([^:\s]+)(?::([\d,]+)-([\d,]+))?$")


def parse_region(region, lengths):
    """'CONTIG' or 'CONTIG:START-END' (1-based, inclusive, commas allowed) -> (contig, start, end)."""
    m = _REGION.match(region.strip())
    example = f"{next(iter(lengths))}:1-100000" if lengths else "CONTIG:1-100000"
    if not m or m.group(1) not in lengths:
        raise ViewerError(f"Region {region!r} not found in the reference; use CONTIG or CONTIG:START-END, e.g. {example}.")
    contig, length = m.group(1), lengths[m.group(1)]
    if m.group(2) is None:
        return contig, 1, length
    start, end = int(m.group(2).replace(",", "")), int(m.group(3).replace(",", ""))
    if not 1 <= start <= end or start > length:
        raise ViewerError(f"Region {region!r} is outside {contig} (1-{length:,}).")
    return contig, start, min(end, length)


# ---- step 2: read a run into page data ----

def _read_manifest(run_dir):
    path = os.path.join(run_dir, INPUTS_DIR, MANIFEST)
    if not os.path.exists(path):
        raise ViewerError(
            f"No viewer inputs in {os.path.join(run_dir, INPUTS_DIR)}. The run's cleanup deletes the "
            "reference FASTA and repeat masks; rerun `coral run_single` with --viewer to keep them.")
    with open(path) as fh:
        return json.load(fh)


def _load_calls(run_dir, taxon, other, reference, windows):
    """{contig: sorted [(position, mutation)]} for `taxon`'s branch inside `windows`."""
    path = os.path.join(run_dir, "Mutations", f"{taxon}__{other}__{reference}__mutations.csv.gz")
    if not os.path.exists(path):
        raise ViewerError(f"Missing {path}; is {run_dir} a finished run_single run for these taxa?")
    calls = {c: [] for c in windows}
    with gzip.open(path, "rt") as fh:
        for row in csv.DictReader(fh):
            window = windows.get(row["chromosome"])
            if window is None:
                continue
            pos = int(row["position"])
            if window[0] <= pos <= window[1]:
                calls[row["chromosome"]].append((pos, row["mutation"]))
    for rows in calls.values():
        rows.sort()
    return calls


def sister_masked_flags(bam_path, sister_mask, contig, positions):
    """For each sorted 1-based call `position` on `contig`: True when every primary read
    covering it came from sequence inside `sister_mask`. Reads are CORAL fragments named
    `<chrom>_<start>_<end>` (1-based, forward strand of the genome they were cut from), so
    a read's query offset maps straight back to a coordinate in that genome."""
    positions = np.asarray(positions, dtype=np.int64)
    n = len(positions)
    if n == 0 or not sister_mask:
        return np.zeros(n, bool)
    reads, masked = np.zeros(n, np.int64), np.zeros(n, np.int64)
    with pysam.AlignmentFile(bam_path) as bam:
        for read in bam.fetch(contig, int(positions[0]) - 1, int(positions[-1])):
            if read.is_unmapped or read.is_secondary:
                continue
            lo = positions.searchsorted(read.reference_start + 1, "left")
            hi = positions.searchsorted(read.reference_end, "right")
            if lo >= hi:
                continue
            try:
                chrom, start, end = read.query_name.rsplit("_", 2)
                start, end = int(start), int(end)
            except ValueError:
                continue
            pairs = np.array(read.get_aligned_pairs(matches_only=True), dtype=np.int64)
            if not len(pairs):
                continue
            ref_pos = pairs[:, 1] + 1
            idx = np.minimum(positions.searchsorted(ref_pos), n - 1)
            hit = positions[idx] == ref_pos
            if not hit.any():
                continue
            cig = read.cigartuples
            hard = cig[0][1] if cig and cig[0][0] == 5 else 0
            offset = pairs[hit, 0] + hard
            length = end - start + 1
            genome_pos = start + (length - 1 - offset if read.is_reverse else offset)
            np.add.at(reads, idx[hit], 1)
            np.add.at(masked, idx[hit], sister_mask.contains_array(chrom, genome_pos).astype(np.int64))
    return (reads > 0) & (masked == reads)


def _coral_version():
    try:
        from importlib.metadata import version
        return version("coral")
    except Exception:
        return None


def build_view_data(run_dir, compare=None, region=None, max_genome_bp=DEFAULT_MAX_GENOME_BP):
    """Everything the page embeds, as a JSON-ready dict."""
    manifest = _read_manifest(run_dir)
    inputs = os.path.join(run_dir, INPUTS_DIR)
    reference, taxa = manifest["reference"], manifest["taxa"]
    params = manifest.get("parameters", {})
    masks = manifest.get("masks", {})
    if len(taxa) != 2:
        raise ViewerError(f"The viewer expects a run_single run with two ingroup taxa; the manifest lists {len(taxa)}.")
    if compare and not params.get("repeat_mask"):
        raise ViewerError("--compare shows the calls a repeat mask removed, but this run was made without --repeat-mask.")

    fasta = pysam.FastaFile(os.path.join(inputs, manifest["fasta"]))
    try:
        lengths = dict(zip(fasta.references, fasta.lengths))
        if region:
            contig, start, end = parse_region(region, lengths)
            windows = {contig: (start, end)}
        else:
            windows = {c: (1, n) for c, n in lengths.items()}
        total = sum(e - s + 1 for s, e in windows.values())
        if total > max_genome_bp:
            largest = max(lengths, key=lengths.get)
            raise ViewerError(
                f"{total:,} bp is more than the viewer embeds (--max-genome-bp {max_genome_bp:,}). "
                f"Build one region instead, e.g. --region {largest}:1-1000000, or raise --max-genome-bp.")

        def load_mask(name):
            return RepeatMask.from_bed(os.path.join(inputs, masks[name])) if name in masks else RepeatMask()

        ref_mask = load_mask(reference)
        has_sister = all(t in masks for t in taxa)
        sister_masks = {t: load_mask(t) for t in taxa}

        kept = {t: _load_calls(run_dir, t, taxa[1 - i], reference, windows) for i, t in enumerate(taxa)}
        removed = {t: {c: [] for c in windows} for t in taxa}
        if compare:
            for i, t in enumerate(taxa):
                for contig, rows in _load_calls(compare, t, taxa[1 - i], reference, windows).items():
                    have = {p for p, _ in kept[t][contig]}
                    removed[t][contig] = [(p, m) for p, m in rows if p not in have]

        contexts, index = [], {}

        def pack_calls(rows, flags):
            pos = np.array([p for p, _ in rows], dtype=np.int64)
            ctx = []
            for _, label in rows:
                if label not in index:
                    index[label] = len(contexts)
                    contexts.append(label)
                ctx.append(index[label])
            return {"d": np.diff(np.r_[0, pos]).tolist() if len(pos) else [],
                    "c": ctx, "f": np.asarray(flags, dtype=np.int64).tolist()}

        contigs = []
        for contig, (start, end) in windows.items():
            packed, n_runs = pack_sequence(fasta.fetch(contig, start - 1, end))
            ms, me = ref_mask.intervals(contig, start, end)
            entry = {
                "name": contig, "length": lengths[contig], "start": start, "end": end,
                "mask": np.column_stack([np.diff(np.r_[0, ms]), me - ms + 1]).ravel().tolist() if len(ms) else [],
                "seq": packed,
                "n": [v + start if k % 2 == 0 else v for k, v in enumerate(n_runs)],   # absolute 1-based starts
                "kept": {}, "removed": {},
            }
            for t in taxa:
                rows = kept[t][contig]
                pos = np.array([p for p, _ in rows], dtype=np.int64)
                flags = np.zeros(len(pos), bool)
                if has_sister and len(pos):
                    bam = os.path.join(run_dir, "BAMs", f"{t}_to_{reference}.bam")
                    if not os.path.exists(bam):
                        raise ViewerError(f"Missing {bam}, needed to trace calls back to {t}'s mask.")
                    flags = sister_masked_flags(bam, sister_masks[t], contig, pos)
                entry["kept"][t] = pack_calls(rows, flags)

                rrows = removed[t][contig]
                rpos = np.array([p for p, _ in rrows], dtype=np.int64)
                near_ref_mask = (ref_mask.contains_array(contig, rpos - 1) | ref_mask.contains_array(contig, rpos)
                                 | ref_mask.contains_array(contig, rpos + 1))
                entry["removed"][t] = pack_calls(rrows, ~near_ref_mask)     # 1 = removed by the sister mask
            contigs.append(entry)
    finally:
        fasta.close()

    meta = {
        "run_id": os.path.basename(os.path.normpath(run_dir)),
        "reference": reference,
        "taxa": taxa,
        "accessions": manifest.get("accessions", {}),
        "parameters": params,
        "masker": params.get("repeat_masker") if params.get("repeat_mask") else None,
        "has_mask": reference in masks,
        "has_sister_masks": has_sister,
        "compare": os.path.basename(os.path.normpath(compare)) if compare else None,
        "region": region,
        "coral_version": _coral_version(),
    }
    return {"meta": meta, "branches": taxa, "contexts": contexts, "contigs": contigs}


def write_html(data, output):
    """Embed `data` in the packaged template and write it to `output` atomically."""
    template = resources.files("coral").joinpath(TEMPLATE).read_text(encoding="utf-8")
    if "__DATA__" not in template:
        raise ViewerError(f"{TEMPLATE} has no __DATA__ placeholder")
    payload = json.dumps(data, separators=(",", ":")).replace("</", "<\\/")
    os.makedirs(os.path.dirname(os.path.abspath(output)), exist_ok=True)
    tmp = output + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        fh.write(template.replace("__DATA__", payload, 1))
    os.replace(tmp, output)
    return output


def build_viewer(run_dir, output=None, compare=None, region=None,
                 max_genome_bp=DEFAULT_MAX_GENOME_BP, verbose=True):
    """Build the page for `run_dir`; returns its path. Default output:
    <run_dir>/Viewer/<run_id>[__<region>].html"""
    log(f"Building viewer for {run_dir}" + (f" (compared with {compare})" if compare else ""), verbose)
    data = build_view_data(run_dir, compare=compare, region=region, max_genome_bp=max_genome_bp)
    if output is None:
        tag = "__" + re.sub(r"[^\w.-]+", "_", region.replace(",", "")) if region else ""
        output = os.path.join(run_dir, VIEWER_DIR, f"{data['meta']['run_id']}{tag}.html")
    write_html(data, output)
    log(f"Viewer written to: {output} ({os.path.getsize(output) / 1e6:.1f} MB)", verbose)
    return output
