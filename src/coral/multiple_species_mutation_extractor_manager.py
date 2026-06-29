import csv
import os
import gzip
import json
from collections import defaultdict
import pandas as pd
from .multiple_species_utils import annotate_tree_with_indices, save_annotated_tree, collapse_mutations, filter_mutations_dict
from .plot_utils import MutationSpectraPlotter
from .utils import log


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
        
        #with open(os.path.join(self.output_dir, "species_mapping.json"), 'w') as f:
        #    json.dump(self.mapping, f, indent=2)
        
        os.makedirs(self.output_dir, exist_ok=True)
        self.plots_dir = os.path.join(self.output_dir, "Plots")
        self.csv_dir = os.path.join(self.output_dir, "CSVs")

    def _all_same(self, seq):
        return len(seq) > 0 and all(ch == seq[0] for ch in seq)

    def _quality_check(self, fields):
        sample_fields = fields[3:]
        return (
            sample_fields
            and all('*' not in field for field in sample_fields)
            and all(self._all_same(field.translate(str.maketrans('', '', '^$[]'))) for field in sample_fields)
        )
    
    def _parse_line(self, line):
        parts = line.strip().split('\t')
        if len(parts) < self.n_species * 3:
            return None
        chrom, pos, ref_base = parts[:3]
        base_calls = parts[4::3]
        normalized = [base[0] if base and base[0] not in {',', '.'} else ref_base for base in base_calls]
        return [chrom, pos, ref_base] + normalized

    def _detect_mutations(self, buffer):
        triplets = [fields[3:] for fields in buffer]
        prev_bases, curr_bases, next_bases = triplets
        if self._all_same(prev_bases) and self._all_same(next_bases) and len(set(curr_bases)) > 1:
            return [
                buffer[1][0],  # chrom
                buffer[1][1],  # pos
                prev_bases[0].upper(),
                next_bases[0].upper(),
                buffer[1][2].upper()
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

    def extract(self):
        csv_path = os.path.join(self.output_dir, "matching_bases.csv.gz")
        #####header = ["chromosome", "position", "left", "right"] + [f"taxa{i}" for i in range(self.n_species)]
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
                        if all(qc_flags):
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

