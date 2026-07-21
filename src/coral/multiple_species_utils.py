from collections import defaultdict
import gzip
import random
import re
import pandas as pd
from io import StringIO
from ete3 import Tree
import sys
import os
import json
from .utils import log

def parse_species_accession_from_newick(newick_str):
    tree = Tree(newick_str, format=1)
    species_accession_dict = {}
    for leaf in tree.iter_leaves():
        if "|" in leaf.name:
            species, accession = leaf.name.split("|", 1)
            species_accession_dict[species] = accession
        else:
            print(f"Leaf name '{leaf.name}' does not contain a '|' separator.")
            sys.exit(1)

    # Determine outgroup: any direct child of root with a single leaf
    outgroup = None
    for child in tree.children:
        leaves = child.get_leaves()
        if len(leaves) == 1:
            outgroup = leaves[0].name.split("|", 1)[0]
            break
    
    if outgroup is None:
        print("Could not determine a single outgroup from the Newick tree. Please ensure the tree is rooted and has a single outgroup.")
        sys.exit(1)

    return species_accession_dict, outgroup


def annotate_tree_with_indices(newick_str, outgroup_name, file_path=None, verbose=True):
    tree = Tree(newick_str, format=1)

    # Normalize leaf names
    for leaf in tree.iter_leaves():
        if '|' in leaf.name:
            leaf.name = leaf.name.split('|', 1)[0]

    # Outgroup first, then rest in tree leaf order (must match pileup BAM order)
    terminals = tree.get_leaves()
    sorted_terminals = [t for t in terminals if t.name == outgroup_name] + \
                       [t for t in terminals if t.name != outgroup_name]

    terminal_mapping = {}
    for idx, node in enumerate(sorted_terminals):
        node.add_feature("index", idx)
        node.add_feature("custom_name", node.name)
        terminal_mapping[idx] = node.name
        terminal_mapping[node.name] = idx

    next_internal_idx = len(sorted_terminals)
    for node in tree.traverse("postorder"):
        if not node.is_leaf():
            node.add_feature("index", next_internal_idx)
            node.add_feature("custom_name", f"Node{next_internal_idx}")
            next_internal_idx += 1
            

    if file_path is not None:
        original_names = {}
        for node in tree.traverse():
            original_names[node] = node.name
            node.name = getattr(node, "custom_name", node.name)

        annotated_tree_path = f"{os.path.splitext(file_path)[0]}_annotated.nwk"
        tree.write(format=1, outfile=annotated_tree_path)

        for node in tree.traverse():
            node.name = original_names[node]

        mapping_path = f"{os.path.splitext(file_path)[0]}_mapping2.json" ### edited
        with open(mapping_path, "w") as f:
            json.dump(terminal_mapping, f, indent=2)

        log(f"Annotated tree saved to {annotated_tree_path}", verbose)
        log(f"Terminal mapping saved to {mapping_path}", verbose)

    sorted_terminal_names = [node.name for node in sorted_terminals]
    return tree, terminal_mapping, sorted_terminal_names


def annotate_list_with_indices(species_list, outgroup_name, file_path=None, verbose=True):

    species_list = [species[0] for species in species_list]
    
    # Outgroup first, then rest in input order (must match pileup BAM order)
    sorted_terminals = [outgroup_name] + [s for s in species_list if s != outgroup_name]

    terminal_mapping = {}
    for idx, node in enumerate(sorted_terminals):
        terminal_mapping[idx] = node
        terminal_mapping[node] = idx

    if file_path is not None:
        mapping_path = f"{os.path.splitext(file_path)[0]}_mapping2.json" ### edited
        with open(mapping_path, "w") as f:
            json.dump(terminal_mapping, f, indent=2)
        
        log(f"Terminal mapping saved to {mapping_path}", verbose)

    return sorted_terminals, terminal_mapping

# Add functions to create mutation spectra from PHYLIP output files.
# and to collapse complementary mutations into a single representation.

def parse_phylip_edges(outfile_path):
    """Edges from the 'between / and / length' table of a dnapars outfile."""
    text = open(outfile_path).read()
    m = re.search(r"between\s+and\s+length\s*\n(.*?)\n\s*\n", text, re.S)
    if m is None:
        raise ValueError(f"No 'between/and/length' table in {outfile_path} "
                         f"(did you pass the .outtree instead of the .outfile?)")
    edges = []
    for line in m.group(1).strip().splitlines():
        p = line.split()
        # keep only real endpoints: an integer (interior) or 'taxaN' (tip)
        if len(p) >= 2 and re.fullmatch(r"\d+|taxa\d+", p[0]) and re.fullmatch(r"\d+|taxa\d+", p[1]):
            edges.append((p[0], p[1]))
    return edges

def phylip_interior_clades(edges, outgroup_label):
    """{interior_number: frozenset(descendant tip labels)} when rooted on the outgroup tip."""
    adj = defaultdict(list)
    for a, b in edges:
        adj[a].append(b); adj[b].append(a)
    clades = {}
    def dfs(node, parent):
        s = {node} if node.startswith("taxa") else set()
        for nb in adj[node]:
            if nb != parent:
                s |= dfs(nb, node)
        clades[node] = frozenset(s)
        return s
    dfs(outgroup_label, None)
    return {n: c for n, c in clades.items() if not n.startswith("taxa")}


def tree_from_phylip_outtree(outtree_path, terminal_mapping, outgroup_name):
    """
    Turn PHYLIP's inferred tree (leaves labeled 'taxa0'..'taxaN', UNROOTED)
    into a tree annotated for Fitch, WITHOUT re-sorting indices — so it stays
    aligned with the columns already in matching_bases.csv.gz.

    terminal_mapping must be the in-memory mapping (int idx -> name AND
    name -> int idx), i.e. self.terminal_mapping, not the JSON-reloaded one.
    """
    tree = Tree(outtree_path, format=1)

    # 1) taxaN  ->  species name
    for leaf in tree.iter_leaves():
        idx = int(leaf.name.replace("taxa", ""))
        leaf.name = terminal_mapping[idx]

    # 2) PHYLIP tree is unrooted -> root on the outgroup so Fitch has polarity
    tree.set_outgroup(tree & outgroup_name)

    # 3) annotate using the EXISTING mapping (do NOT call annotate_tree_with_indices() here,
    #    because it sorts alphabetically and would desync leaf indices from the CSV columns)
    for leaf in tree.iter_leaves():
        leaf.add_feature("index", terminal_mapping[leaf.name])
        leaf.add_feature("custom_name", leaf.name)

    """
    next_idx = sum(1 for _ in tree.iter_leaves())
    for node in tree.traverse("postorder"):
        if not node.is_leaf():
            node.add_feature("index", next_idx)
            node.add_feature("custom_name", f"Node({next_idx})")
            next_idx += 1
    """
    outfile_path = outtree_path.replace(".outtree", ".outfile")   # the table lives in .outfile
    edges = parse_phylip_edges(outfile_path)
    clades = phylip_interior_clades(edges, f"taxa{terminal_mapping[outgroup_name]}")
    next_idx = sum(1 for _ in tree.iter_leaves())
    for node in tree.traverse("postorder"):
        if not node.is_leaf():
            clade = frozenset(f"taxa{terminal_mapping[l.name]}" for l in node.iter_leaves())
            match = next((n for n, c in clades.items() if c == clade), None)
            node.add_feature("index", next_idx); next_idx += 1        # keep a unique index
            node.add_feature("custom_name",
                             match if match is not None
                             else ("ROOT" if node.is_root() else f"Node({next_idx-1})"))

    return tree

def save_annotated_tree(tree, path):
    original_names = {}
    for node in tree.traverse():
        original_names[node] = node.name
        node.name = getattr(node, "custom_name", node.name)

    tree.write(format=1, outfile=path)

    for node in tree.traverse():
        node.name = original_names[node]

complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

def get_complement(mutation):
    comp = [complement[nuc] if nuc in complement else nuc for nuc in mutation]
    comp[0], comp[-1] = comp[-1], comp[0]
    return ''.join(comp)

def collapse_mutations(mutation_dict):
    collapsed = defaultdict(int)
    for mutation, count in mutation_dict.items():
        if mutation[2] in {'A', 'G'}:
            collapsed[get_complement(mutation)] += int(count)
        else:
            collapsed[mutation] += int(count)
    return collapsed


def load_random_rows(file_path, max_rows=1000000, seed=42, verbose=True):
    random.seed(seed)
    
    # Count total rows (excluding header)
    with gzip.open(file_path, 'rt') as f:
        header = f.readline()
        total_rows = sum(1 for _ in f)
    
    log(f"File has {total_rows} rows (excluding header).", verbose)

    if total_rows <= max_rows:
        log("Loading full gzipped file.", verbose)
        return pd.read_csv(file_path, index_col=0, compression='gzip').astype(str)
    
    sampled_indices = set(random.sample(range(total_rows), max_rows))

    with gzip.open(file_path, 'rt') as f:
        header = f.readline()
        sampled_lines = [line for i, line in enumerate(f) if i in sampled_indices]

    return pd.read_csv(StringIO(header + ''.join(sampled_lines)), index_col=0).astype(str)

mutation_pattern = re.compile(r"^[ACGT]\[[ACGT]>[ACGT]\][ACGT]$")

def filter_mutations_dict(d):
    return {k: v for k, v in d.items() if mutation_pattern.match(k)}




valid_bases = {'A', 'C', 'G', 'T'}

def collapse_triplets(triplet_dict):
    """Pyrimidine-centric fold: reverse-complement a triplet whose CENTER base
    (triplet[1]) is a purine (G/A). Mirror of MutationNormalizer.collapse_triplets
    so collapsed triplet keys line up with collapsed mutation keys."""
    collapsed = defaultdict(int)
    for triplet, count in triplet_dict.items():
        if triplet[1] in {'G', 'A'}:
            rc = ''.join(complement.get(nuc, nuc) for nuc in reversed(triplet))
            collapsed[rc] += int(count)
        else:
            collapsed[triplet] += int(count)
    return collapsed

def filter_triplets_dict(d):
    """Drop triplets containing 'N' or any non-ACGT char. Mirror of
    MutationNormalizer.filter_triplets_dict."""
    return {k: v for k, v in d.items() if 'N' not in k and all(x in valid_bases for x in k)}

def normalize_by_triplets(mutations, triplets):
    """Divide each collapsed mutation count by its trinucleotide context count
    ('X[R>A]Y' -> key 'XRY' = k[0]+k[2]+k[-1]); 0 when context unseen. Mirror of
    MutationNormalizer.normalize_by_triplets."""
    out = {}
    for k, v in mutations.items():
        denom = triplets.get(f"{k[0]}{k[2]}{k[-1]}", 0)
        out[k] = v / denom if denom > 0 else 0
    return out

def scale_counts(d, target_sum=10000):
    """Rescale a count dict to sum to target_sum (0 when empty). Mirror of
    MutationNormalizer.scale_counts."""
    total = sum(d.values())
    if total == 0:
        return {k: 0 for k in d}
    return {k: round(v / total * target_sum) for k, v in d.items()}