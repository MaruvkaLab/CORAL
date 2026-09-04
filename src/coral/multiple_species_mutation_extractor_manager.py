import csv
import json
import os
import gzip
from collections import defaultdict
import pandas as pd
from .multiple_species_utils import annotate_tree_with_indices, save_annotated_tree
from .mutation_extractor_manager import MutationNormalizer
from .plot_utils import MutationSpectraPlotter
from .utils import log
import re

MIN_DEPTH = 1  # Minimum depth threshold for quality check

# Matches read-start (^ + mapQ char), read-end ($), or an indel marker (+N / -N).
# Digits are captured so the N following indel bases can be skipped in one pass.
_CLEAN_RE = re.compile(r'\^.|\$|[+-](\d*)')


def clean_bases(s):
    """Strip read-start/end markers and indel notation from an mpileup bases field.

    Fast path: every removable marker (^x, $, +N…, -N…) contains at least one of
    ^, $, +, or -. When none are present—the overwhelmingly common case for
    CORAL's shallow pileups—the regex would return s unchanged, so we skip it.
    This is byte-for-byte identical to the regex path while avoiding the regex
    call and list construction.
    """
    if '^' not in s and '$' not in s and '+' not in s and '-' not in s:
        return s
    out = []
    pos = 0
    for m in _CLEAN_RE.finditer(s):
        if m.start() < pos:
            continue  # falls inside a just-skipped indel run; can't happen in valid pileup syntax
        out.append(s[pos:m.start()])
        digits = m.group(1)
        if digits is not None:
            pos = m.end() + (int(digits) if digits else 0)
        else:
            pos = m.end()
    out.append(s[pos:])
    return ''.join(out)


class ParsedLine(list):
    """Parsed pileup row: [chrom, pos, ref, (depth, bases), ...].

    Subclasses list to preserve native list-indexing performance.
    clean_sample(i) caches cleaned bases per sample, avoiding repeated
    clean_bases() calls when a row is reused across sliding windows.
    """
    __slots__ = ('_clean',)

    def __init__(self, chrom, pos, ref_base, samples):
        super().__init__([chrom, pos, ref_base] + samples)
        self._clean = [None] * len(samples)

    def clean_sample(self, sample_idx):
        cached = self._clean[sample_idx]
        if cached is None:
            _, bases = self[3 + sample_idx]
            cached = clean_bases(bases)
            self._clean[sample_idx] = cached
        return cached


class MultipleSpeciesMutationExtractor:
    def __init__(self, pileup_file, output_dir, n_species, tree=None, species_list=None, mapping=None, no_cache=False, verbose=False,
                 ref_fasta=None, bams=None, fai_path=None, cores=1, scan_jobs=None, max_memory_mb=None,
                 fitch_jobs=None, parallel_fitch=False):
        self.pileup_file = pileup_file
        self.output_dir = output_dir
        self.n_species = n_species
        self.tree = tree
        self.species_list = species_list
        self.mapping = mapping
        self.no_cache = no_cache
        self.verbose = verbose

        # Optional; omitting these reproduces the serial behaviour exactly.
        self.ref_fasta = ref_fasta
        self.bams = list(bams) if bams else None
        self.fai_path = fai_path
        self.cores = cores or 1
        self.scan_jobs = scan_jobs
        self.max_memory_mb = max_memory_mb
        self.fitch_jobs = fitch_jobs
        self.parallel_fitch = parallel_fitch
        if self.tree is None and self.species_list is None:
            raise ValueError("Either newick_tree or species_list must be provided.")
        if self.mapping is None:
            raise ValueError("Dictionary mapping taxa names must be provided.")

        os.makedirs(self.output_dir, exist_ok=True)
        self.plots_dir = os.path.join(self.output_dir, "Plots")
        self.csv_dir = os.path.join(self.output_dir, "Mutations")
        self.triplet_counts = None  # trinucleotide opportunities, from the scan

    def _all_same(self, seq):
        # count() runs its comparison loop in C; a generator + all() pays
        # per-element Python-level overhead for what is otherwise a tight scan.
        return len(seq) > 0 and seq.count(seq[0]) == len(seq)

    def _parse_line(self, line):
        parts = line.strip().split('\t')
        if len(parts) < self.n_species * 3:
            return None
        chrom, pos, ref_base = parts[:3]
        depths = parts[3::3]
        base_calls = parts[4::3]
        return ParsedLine(chrom, pos, ref_base, list(zip(depths, base_calls)))

    def _quality_check(self, fields):
        if not fields:
            return False
        samples = fields[3:]
        for i, (depth, bases) in enumerate(samples):
            if '*' in bases: # deletions
                return False
            if '+' in bases: # insertions
                return False
            if int(depth) < MIN_DEPTH:
                return False
            cleaned = fields.clean_sample(i).replace(',', '.').lower()
            if not self._all_same(cleaned):
                return False
        return True

    def _detect_site(self, buffer):
        """One look at the 3-line window -> (mutation_row, triplet_context).

        Flanks must match their own ref base. The centre is then one of:
        ingroup varies -> row + counted context; all at ref -> counted context
        only; no sample at ref -> neither, as a single outgroup cannot polarise
        it (same exclusion as the two-taxa detect_mutation_triplet).
        """
        def normalize(fields, ref_base):
            rb = ref_base.upper()
            out = []
            for i in range(len(fields) - 3):
                cleaned = fields.clean_sample(i)
                c = cleaned[0].upper() if cleaned else ''
                out.append(rb if (not c or c in {',', '.'}) else c)
            return out

        def matches_ref(fields):
            rb = fields[2].upper()
            return all(b == rb for b in normalize(fields, rb))

        if not (matches_ref(buffer[0]) and matches_ref(buffer[2])):
            return None, None

        ref_upper = buffer[1][2].upper()
        curr_bases = normalize(buffer[1], buffer[1][2])   # already uppercased
        distinct = set(curr_bases)

        if ref_upper not in distinct:
            return None, None
        triplet_context = buffer[0][2].upper() + ref_upper + buffer[2][2].upper()

        if len(distinct) > 1:   # variable across the ingroup -> CSV row
            mutation_row = [
                buffer[1][0],           # chrom
                buffer[1][1],           # pos
                buffer[0][2].upper(),   # left  = prev line's ref base
                buffer[2][2].upper(),   # right = next line's ref base
                ref_upper,              # = taxa0, the outgroup (it has no BAM)
            ] + [b.upper() for b in curr_bases]
            return mutation_row, triplet_context

        return None, triplet_context   # fully invariant

    @property
    def parallel_scan(self):
        """True when the scan can run per chromosome: it needs the indexed BAMs
        and the reference, which the serial path does not."""
        return bool(self.ref_fasta and self.bams and self.fai_path)


    def _recursive_state_check(self, node, row):
        if node.is_leaf():
            node.add_feature("state", {row[f"taxa{self.mapping[node.name]}"]})
            return node.state
        # Handle nodes with any number of children (supporting multifurcating trees)
        child_states = [self._recursive_state_check(child, row) for child in node.children]
        # Intersect all child states if any intersection exists, otherwise union
        node_state = child_states[0]
        for child_state in child_states[1:]:
            intersect = node_state & child_state
            node_state = intersect if intersect else node_state | child_state
        node.add_feature("state", node_state)
        return node_state

    def _recursive_fitch(self, node, parent_state, row, mutation_dict, ambiguous_count):
        next_state = parent_state
        # Ensure node.state exists (should be set by _recursive_state_check, but add safety check)
        if not hasattr(node, 'state'):
            raise RuntimeError(f"Node {node.name} missing state attribute. Tree may not have been properly initialized.")
        if parent_state not in node.state:
            if len(node.state) > 1:
                return mutation_dict, ambiguous_count + 1
            next_state = list(node.state)[0]
            parent_name = node.up.custom_name if node.up else "ROOT"
            branch_key = f"{parent_name}→{node.custom_name}"
            mutation = f"{row['left']}[{parent_state}>{next_state}]{row['right']}"
            mutation_dict.setdefault(branch_key, []).append((row['chromosome'], row['position'], mutation))
        if not node.is_leaf():
            for child in node.children:
                mutation_dict, ambiguous_count = self._recursive_fitch(child, next_state, row, mutation_dict, ambiguous_count)
        return mutation_dict, ambiguous_count

    def _fitch(self, tree_root, row, mutation_dict):
        root_state = self._recursive_state_check(tree_root, row)
        if len(root_state) == 1:
            return self._recursive_fitch(tree_root, list(root_state)[0], row, mutation_dict, 0)
        return mutation_dict, 1

    def _consecutive(self, *lines):
        chrom = lines[0][0]                       # index 0 = chrom
        positions = [int(f[1]) for f in lines]    # index 1 = pos
        return (all(f[0] == chrom for f in lines)
                and all(positions[i] + 1 == positions[i + 1] for i in range(len(positions) - 1)))

    def extract(self):
        csv_path = os.path.join(self.output_dir, "matching_bases.csv.gz")
        header = ["chromosome", "position", "left", "right"] + [f"taxa{k}" for k in self.mapping if isinstance(k, int)]

        triplets_path = os.path.join(self.output_dir, "triplets.json")

        # Both scan outputs are required: the scan clears triplets.json before it
        # starts, so a lone CSV means the scan did not finish.
        if os.path.exists(csv_path) and os.path.exists(triplets_path) and not self.no_cache:
            log(f'Using cached matching positions from csv at {csv_path}', self.verbose)
            with open(triplets_path) as tf:
                self.triplet_counts = json.load(tf)
        elif self.parallel_scan:
            from .parallel import scan_pileup_parallel   # keeps mp off the serial path
            self.triplet_counts = scan_pileup_parallel(self, csv_path, triplets_path, header)
        else:
            self.triplet_counts = self._scan_pileup_serial(csv_path, triplets_path, header)

        if self.tree:
            if self.parallel_fitch:
                from .parallel import run_fitch
                mutation_dict, ambiguous_counter = run_fitch(
                    csv_path, self.tree, self.mapping,
                    cores=self.cores,
                    fitch_jobs=self.fitch_jobs,
                    max_memory_mb=self.max_memory_mb,
                    log_fn=(lambda msg: log(msg, self.verbose)),
                )
            else:
                mutation_dict, ambiguous_counter = self._fitch_serial(csv_path)
            self._save_results(mutation_dict)
            log(f"Total ambiguous mutations: {ambiguous_counter}", self.verbose)

    def _scan_pileup_serial(self, csv_path, triplets_path, header):
        """Scan the whole-genome pileup with a 3-line sliding window.

        Writes both scan outputs and returns the triplet counts, matching
        parallel.scan_pileup_parallel. Written to a temporary name and renamed
        so failed jobs don't leave truncated files.
        """
        triplet_counts = defaultdict(int)
        tmp_path = csv_path + ".tmp"
        # One generation with the CSV: clear stale counts so a crash between the
        # two writes can't pair a fresh CSV with old triplets.
        if os.path.exists(triplets_path):
            os.remove(triplets_path)
        try:
            with gzip.open(tmp_path, 'wt', newline='') as outfile:
                writer = csv.writer(outfile)
                writer.writerow(header)

                with gzip.open(self.pileup_file, 'rt') as infile:
                    buffer = [None, self._parse_line(infile.readline()), self._parse_line(infile.readline())]
                    qc_flags = [False, self._quality_check(buffer[1]), self._quality_check(buffer[2])]

                    for line in infile:
                        buffer = [buffer[1], buffer[2], self._parse_line(line)]
                        qc_flags = [qc_flags[1], qc_flags[2], self._quality_check(buffer[2])]
                        if all(qc_flags) and self._consecutive(*buffer):
                            result, triplet = self._detect_site(buffer)
                            if triplet is not None:
                                triplet_counts[triplet] += 1
                            if result is not None:
                                writer.writerow(result)
            os.replace(tmp_path, csv_path)
        except BaseException:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)
            raise

        triplet_counts = dict(triplet_counts)
        # Renamed like the CSV, so a kill mid-write can't leave truncated JSON.
        tmp_triplets = triplets_path + ".tmp"
        with open(tmp_triplets, "w") as tf:
            json.dump(triplet_counts, tf, indent=2)
        os.replace(tmp_triplets, triplets_path)
        log(f"Saved {len(triplet_counts)} triplet contexts to {triplets_path}", self.verbose)
        return triplet_counts

    def _fitch_serial(self, csv_path):
        """Fitch reconstruction, one row at a time over the scanned CSV."""
        mutation_dict = defaultdict(list)
        ambiguous_counter = 0

        for chunk in pd.read_csv(csv_path, chunksize=1000):
            for _, row in chunk.iterrows():
                mutation_dict, ambiguous = self._fitch(self.tree.copy(), row, mutation_dict)
                ambiguous_counter += ambiguous

        return mutation_dict, ambiguous_counter

    def _save_results(self, mutation_dict):
        spectra_plotter = MutationSpectraPlotter()
        os.makedirs(self.plots_dir, exist_ok=True)
        os.makedirs(self.csv_dir, exist_ok=True)
        # run_single's normaliser, so both pipelines fold and scale identically.
        # It owns Tables/ and creates it.
        normalizer = MutationNormalizer(self.output_dir, verbose=self.verbose)
        tables_dir = normalizer.output_dir

        # One shared, branch-agnostic denominator for every branch.
        collapsed_triplets = normalizer.collapse_triplets(
            normalizer.filter_triplets_dict(self.triplet_counts or {}))

        spectra_dict = {}       # per-branch collapsed+filtered raw spectrum
        normalized_scaled = {}  # per-branch spectrum / triplet opportunities, scaled
        scaled_raw = {}         # per-branch raw spectrum, scaled

        for branch_key, mutations in mutation_dict.items():
            df = pd.DataFrame(mutations, columns=["chromosome", "position", "mutation"])
            csv_path = os.path.join(self.csv_dir, f"{branch_key}.csv.gz")
            df.to_csv(csv_path, index=False, header=False, sep="\t", compression="gzip")
            mutation_spectra = normalizer.collapse_mutations(dict(df['mutation'].value_counts()))
            mutation_spectra = normalizer.filter_mutations_dict(mutation_spectra)
            spectra_dict[branch_key] = mutation_spectra
            normalized_scaled[branch_key] = normalizer.scale_counts(
                normalizer.normalize_by_triplets(mutation_spectra, collapsed_triplets))
            scaled_raw[branch_key] = normalizer.scale_counts(mutation_spectra)
            spectra_plot_path = os.path.join(self.plots_dir, f"{branch_key}_spectra.png")
            spectra_plotter.plot_mutations(pd.Series(mutation_spectra), spectra_plot_path, f"Mutation Spectra: {branch_key}")

        spectra_df = pd.DataFrame(spectra_dict)
        spectra_df.to_csv(os.path.join(self.output_dir, "mutation_spectras.tsv"), sep="\t")

        spectra_df.to_csv(os.path.join(tables_dir, "collapsed_mutations.tsv"), sep="\t")
        pd.DataFrame(normalized_scaled).to_csv(os.path.join(tables_dir, "normalized_scaled.tsv"), sep="\t")
        pd.DataFrame(scaled_raw).to_csv(os.path.join(tables_dir, "scaled_raw.tsv"), sep="\t")
        # The denominator is one shared vector, so a single column, not per-branch.
        pd.Series(collapsed_triplets, name="triplets").to_frame().to_csv(
            os.path.join(tables_dir, "triplets.tsv"), sep="\t")
