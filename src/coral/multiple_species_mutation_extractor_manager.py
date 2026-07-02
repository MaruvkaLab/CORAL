import csv
import os
import gzip
from collections import defaultdict
import pandas as pd
from .multiple_species_utils import annotate_tree_with_indices, save_annotated_tree, collapse_mutations, filter_mutations_dict
from .plot_utils import MutationSpectraPlotter
from .utils import log
import re

MIN_DEPTH = 1  # Minimum depth threshold for quality check

# Matches read-start (^ + mapQ char), read-end ($), or an indel marker (+N / -N).
# Digits are captured so the N following indel bases can be skipped in one pass.
_CLEAN_RE = re.compile(r'\^.|\$|[+-](\d*)')


def clean_bases(s):
    """Strip read-start/end markers and indel notation from an mpileup bases field.

    Single regex scan instead of a per-character Python loop: for the common
    case (no indels/read boundaries at this position) finditer finds nothing
    and this degrades to one slice, at roughly C speed.
    """
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
    """A parsed pileup row: [chrom, pos, ref, (depth, bases), (depth, bases), ...].
    Subclasses list (not composes over it) so fields[0]/fields[3:]/etc. stay
    native C-level list indexing — no per-access Python method-call overhead.
    Adds one thing: clean_sample(i) caches clean_bases() per sample so QC and
    mutation detection never clean the same bases string twice. A line can
    appear in up to three sliding-window triplets, so without caching it
    could get re-cleaned 3-4x over."""
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
    def __init__(self, pileup_file, output_dir, n_species, tree=None, species_list=None, mapping=None, no_cache=False, verbose=False):
        self.pileup_file = pileup_file
        self.output_dir = output_dir
        self.n_species = n_species
        self.tree = tree
        self.species_list = species_list
        self.mapping = mapping
        self.no_cache = no_cache
        self.verbose = verbose
        if self.tree is None and self.species_list is None:
            raise ValueError("Either newick_tree or species_list must be provided.")
        if self.mapping is None:
            raise ValueError("Dictionary mapping taxa names must be provided.")

        os.makedirs(self.output_dir, exist_ok=True)
        self.plots_dir = os.path.join(self.output_dir, "Plots")
        self.csv_dir = os.path.join(self.output_dir, "CSVs")

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
    
    """
    def _detect_mutations(self, buffer):
        def normalize(fields, ref_base):
            rb = ref_base.upper()
            out = []
            for i in range(len(fields[3:])):
                cleaned = fields.clean_sample(i)
                c = cleaned[0].upper() if cleaned else ''
                out.append(rb if (not c or c in {',', '.'}) else c)
            return out

        ref_base = buffer[1][2]
        prev_bases = normalize(buffer[0], buffer[0][2])
        curr_bases = normalize(buffer[1], ref_base)
        next_bases = normalize(buffer[2], buffer[2][2])

        if self._all_same(prev_bases) and self._all_same(next_bases) and len(set(curr_bases)) > 1:
            return [
                buffer[1][0], # chrom
                buffer[1][1], # pos
                prev_bases[0].upper(),
                next_bases[0].upper(),
                ref_base.upper()
            ] + [b.upper() for b in curr_bases]
        return None
    """

    def _detect_mutations(self, buffer):
        def normalize(fields, ref_base):
            rb = ref_base.upper()
            out = []
            for i in range(len(fields[3:])):
                cleaned = fields.clean_sample(i)
                c = cleaned[0].upper() if cleaned else ''
                out.append(rb if (not c or c in {',', '.'}) else c)
            return out

        def matches_ref(fields):
            rb = fields[2].upper()
            return all(b == rb for b in normalize(fields, rb))

        ref_base = buffer[1][2]
        curr_bases = normalize(buffer[1], ref_base)

        if matches_ref(buffer[0]) and matches_ref(buffer[2]) and len(set(curr_bases)) > 1:
            return [
                buffer[1][0],  # chrom
                buffer[1][1],  # pos
                buffer[0][2].upper(),  # left context = prev line's ref base
                buffer[2][2].upper(),  # right context = next line's ref base
                ref_base.upper()
            ] + [b.upper() for b in curr_bases]
        return None
    
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

        if os.path.exists(csv_path) and not self.no_cache:
            log(f'Using cached matching positions from csv at {csv_path}', self.verbose)
        else:
            with gzip.open(csv_path, 'wt', newline='') as outfile:
                writer = csv.writer(outfile)
                writer.writerow(header)

                with gzip.open(self.pileup_file, 'rt') as infile:
                    buffer = [None, self._parse_line(infile.readline()), self._parse_line(infile.readline())]
                    qc_flags = [False, self._quality_check(buffer[1]), self._quality_check(buffer[2])]

                    for line in infile:
                        buffer = [buffer[1], buffer[2], self._parse_line(line)]
                        qc_flags = [qc_flags[1], qc_flags[2], self._quality_check(buffer[2])]
                        if all(qc_flags) and self._consecutive(*buffer):
                            result = self._detect_mutations(buffer)
                            if result:
                                writer.writerow(result)

        if self.tree:
            mutation_dict = defaultdict(list)
            ambiguous_counter = 0

            for chunk in pd.read_csv(csv_path, chunksize=1000):
                for _, row in chunk.iterrows():
                    mutation_dict, ambiguous = self._fitch(self.tree.copy(), row, mutation_dict)
                    ambiguous_counter += ambiguous

            self._save_results(mutation_dict)
            log(f"Total ambiguous mutations: {ambiguous_counter}", self.verbose)

    def _save_results(self, mutation_dict):
        spectra_plotter = MutationSpectraPlotter()
        os.makedirs(self.plots_dir, exist_ok=True)
        os.makedirs(self.csv_dir, exist_ok=True)
        spectra_dict = {}

        for branch_key, mutations in mutation_dict.items():
            df = pd.DataFrame(mutations, columns=["chromosome", "position", "mutation"])
            csv_path = os.path.join(self.csv_dir, f"{branch_key}.csv.gz")
            df.to_csv(csv_path, index=False, header=False, sep="\t", compression="gzip")
            mutation_spectra = collapse_mutations(dict(df['mutation'].value_counts()))
            mutation_spectra = filter_mutations_dict(mutation_spectra)
            spectra_dict[branch_key] = mutation_spectra
            spectra_plot_path = os.path.join(self.plots_dir, f"{branch_key}_spectra.png")
            spectra_plotter.plot_mutations(pd.Series(mutation_spectra), spectra_plot_path, f"Mutation Spectra: {branch_key}")

        spectra_df = pd.DataFrame(spectra_dict)
        spectra_df.to_csv(os.path.join(self.output_dir, "mutation_spectras.tsv"), sep="\t")
