from collections import defaultdict
import gzip
import random
import re
import pandas as pd
from io import StringIO
import sys
import os
import json
from .utils import log

def parse_species_accession_from_newick(newick_str):
    from ete3 import Tree
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
    from ete3 import Tree
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
            node.add_feature("custom_name", f"Node({next_internal_idx})")
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


def save_annotated_tree(tree, path):
    original_names = {}
    for node in tree.traverse():
        original_names[node] = node.name
        node.name = getattr(node, "custom_name", node.name)

    tree.write(format=1, outfile=path)

    for node in tree.traverse():
        node.name = original_names[node]

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


def read_phylip_outtree_newick(outtree_path, verbose=True):
    r"""The first newick string in a PHYLIP .outtree. """
    text = open(outtree_path).read()

    start = text.find("(")
    if start == -1:
        raise ValueError(f"No newick tree found in {outtree_path}")
    end = text.find(";", start)
    if end == -1:
        raise ValueError(f"Unterminated newick tree in {outtree_path} "
                         f"(no ';' after the opening '(')")

    if text[end + 1:].strip():
        log(f"{outtree_path} holds more than one tree; using the first.", verbose)

    return "".join(text[start:end + 1].split())


def tree_from_phylip_outtree(outtree_path, terminal_mapping, outgroup_name):
    """PHYLIP's tree (leaves 'taxa0'..'taxaN', unrooted) annotated for Fitch.

    Indices come from `terminal_mapping`, so leaves stay
    aligned with the taxaK columns in matching_bases.csv.gz.
    """
    from ete3 import Tree

    newick = read_phylip_outtree_newick(outtree_path)
    try:
        tree = Tree(newick, format=1)
    except Exception as exc:
        raise ValueError(
            f"Could not parse the PHYLIP tree in {outtree_path}: {exc}\n"
            f"newick read: {newick[:200]}"
        ) from exc

    for leaf in tree.iter_leaves():                    # taxaN -> species name
        idx = int(leaf.name.replace("taxa", ""))
        leaf.name = terminal_mapping[idx]

    tree.set_outgroup(tree & outgroup_name)            # unrooted -> polarity for Fitch

    for leaf in tree.iter_leaves():
        leaf.add_feature("index", terminal_mapping[leaf.name])
        leaf.add_feature("custom_name", leaf.name)

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
