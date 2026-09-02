"""Fitch reconstruction without the per-row tree copy, optionally over a pool.

CORAL's multi-species extractor runs Fitch as::

    for chunk in pd.read_csv(csv_path, chunksize=1000):
        for _, row in chunk.iterrows():
            mutation_dict, ambiguous = self._fitch(self.tree.copy(), row, mutation_dict)

``self.tree.copy()`` is an ete3 deep copy paid *per site*, only because
``_recursive_state_check`` mutates the tree by attaching a ``state`` feature to
every node. Two independent changes are made here, both output-preserving:

**(a) Flatten the tree once.** The topology is decomposed a single time into
plain arrays (``children``, ``custom_name``, ``parent_name``, and per leaf its
``taxaK`` column). Fitch then writes its states into a per-row scratch list, so
nothing about the tree is mutated and no copy is needed.

**(b) Parallelise over row blocks.** The parent process keeps reading with
``pd.read_csv(..., chunksize=N)`` -- deliberately, so dtypes (and hence how
``position`` is rendered in ``Mutations/*.csv.gz``) are exactly what the serial
path produced -- and streams the chunks through ``imap``. ``imap`` yields in
input order and chunks are contiguous row blocks, so concatenating each branch's
lists in chunk order reproduces the serial record order exactly.

The recursion below is a transcription of ``_recursive_state_check`` and
``_recursive_fitch``: same intersect-else-union fold, same child order, same
early ``return`` on an ambiguous node (which skips that whole subtree), same
branch-key and mutation-string formatting, same ambiguity accounting.
"""

from __future__ import annotations

import pandas as pd

from .parallel_resources import cap_jobs, get_mp_context

#: Rows per block handed to a worker. Matches the parent's ``chunksize=1000``
#: read granularity by default; larger blocks amortise pickling on wide trees.
DEFAULT_CHUNK_SIZE = 20_000

#: Rough resident memory per Fitch worker (one DataFrame block + its results).
FITCH_MB_PER_JOB = 800


class FlatTree:
    """A Fitch-ready, picklable, immutable view of an annotated CORAL tree.

    Node ids are assigned in the same pre-order the original recursion visits;
    only the *relative order of each node's children* affects the output, and
    that is preserved from ``node.children``.
    """

    __slots__ = ("children", "custom_name", "parent_name", "leaf_col", "root", "taxa_columns")

    def __init__(self, tree, mapping):
        self.children = []
        self.custom_name = []
        self.parent_name = []
        self.leaf_col = []      # index into the row's taxa values, or -1 for internal
        self.root = 0

        # The ``taxaK`` column order of matching_bases.csv.gz, as built by the
        # extractor: the integer keys of the mapping in insertion order.
        self.taxa_columns = [f"taxa{k}" for k in mapping if isinstance(k, int)]
        col_index = {name: i for i, name in enumerate(self.taxa_columns)}

        def add(node):
            idx = len(self.children)
            self.children.append([])
            self.custom_name.append(node.custom_name)
            self.parent_name.append(node.up.custom_name if node.up else "ROOT")
            if node.is_leaf():
                self.leaf_col.append(col_index[f"taxa{mapping[node.name]}"])
            else:
                self.leaf_col.append(-1)
            for child in node.children:
                self.children[idx].append(add(child))
            return idx

        add(tree)

    def __reduce__(self):
        # __slots__ without __dict__: define pickling explicitly for spawn.
        return (_rebuild_flat_tree, (self.children, self.custom_name,
                                     self.parent_name, self.leaf_col,
                                     self.root, self.taxa_columns))


def _rebuild_flat_tree(children, custom_name, parent_name, leaf_col, root, taxa_columns):
    flat = FlatTree.__new__(FlatTree)
    flat.children = children
    flat.custom_name = custom_name
    flat.parent_name = parent_name
    flat.leaf_col = leaf_col
    flat.root = root
    flat.taxa_columns = taxa_columns
    return flat


def _states(flat, taxa_values):
    """Fitch downward pass: the state set of every node for one site.

    Transcribes ``_recursive_state_check``. Sets are never mutated in place, so a
    single-child node legitimately shares its child's set object (as the original
    does when it assigns ``node_state = child_states[0]``).
    """
    states = [None] * len(flat.children)

    def visit(i):
        col = flat.leaf_col[i]
        if col >= 0:
            states[i] = {taxa_values[col]}
            return states[i]
        child_states = [visit(c) for c in flat.children[i]]
        node_state = child_states[0]
        for child_state in child_states[1:]:
            intersect = node_state & child_state
            node_state = intersect if intersect else node_state | child_state
        states[i] = node_state
        return node_state

    visit(flat.root)
    return states


def _descend(flat, states, i, parent_state, chrom, pos, left, right, out):
    """Fitch upward pass for one site; returns this subtree's ambiguity count."""
    next_state = parent_state
    if parent_state not in states[i]:
        if len(states[i]) > 1:
            # Ambiguous: the original returns here, so the subtree is not walked.
            return 1
        next_state = next(iter(states[i]))
        branch_key = f"{flat.parent_name[i]}→{flat.custom_name[i]}"
        mutation = f"{left}[{parent_state}>{next_state}]{right}"
        out.setdefault(branch_key, []).append((chrom, pos, mutation))
    ambiguous = 0
    for child in flat.children[i]:
        ambiguous += _descend(flat, states, child, next_state, chrom, pos, left, right, out)
    return ambiguous


def fitch_block(flat, block):
    """Run Fitch over one row block. Returns ``(mutation_dict, ambiguous_count)``.

    ``block`` is a DataFrame read by pandas, so ``position`` keeps the dtype the
    serial path had and is rendered identically downstream.
    """
    out = {}
    ambiguous = 0
    n_meta = 4  # chromosome, position, left, right
    columns = ["chromosome", "position", "left", "right"] + flat.taxa_columns
    for row in block[columns].itertuples(index=False, name=None):
        chrom, pos, left, right = row[:n_meta]
        states = _states(flat, row[n_meta:])
        root_state = states[flat.root]
        if len(root_state) != 1:
            ambiguous += 1
            continue
        ambiguous += _descend(flat, states, flat.root, next(iter(root_state)),
                              chrom, pos, left, right, out)
    return out, ambiguous


# -- pool plumbing ---------------------------------------------------------------

_WORKER_FLAT = None


def _init_worker(flat):
    global _WORKER_FLAT
    _WORKER_FLAT = flat


def _fitch_block_worker(block):
    return fitch_block(_WORKER_FLAT, block)


def run_fitch(csv_path, tree, mapping, cores=1, fitch_jobs=None, max_memory_mb=None,
              chunk_size=DEFAULT_CHUNK_SIZE, log_fn=None):
    """Flattened Fitch over ``matching_bases.csv.gz``, in one process or a pool.

    Returns ``(mutation_dict, ambiguous_count)`` ready for ``_save_results``.
    With ``fitch_jobs=1`` this is the flattened single-process implementation
    (still much faster than the per-row ``tree.copy()``); the untouched original
    recursion remains available through the serial pipeline path.
    """
    flat = FlatTree(tree, mapping)
    # Task count is unknown until the file is read; cap on cores and let the pool
    # size be the ceiling rather than the exact block count.
    n_workers = cap_jobs(fitch_jobs, cores, cores,
                         mb_per_job=FITCH_MB_PER_JOB, max_memory_mb=max_memory_mb)

    blocks = pd.read_csv(csv_path, chunksize=chunk_size)
    mutation_dict = {}
    ambiguous_total = 0

    def absorb(result):
        nonlocal ambiguous_total
        out, ambiguous = result
        ambiguous_total += ambiguous
        for branch_key, records in out.items():
            # imap preserves input order and blocks are contiguous, so appending
            # here reproduces the serial per-branch record order.
            mutation_dict.setdefault(branch_key, []).extend(records)

    if n_workers <= 1:
        if log_fn:
            log_fn("Running Fitch reconstruction in a single process (flattened tree)...")
        for block in blocks:
            absorb(fitch_block(flat, block))
    else:
        if log_fn:
            log_fn(f"Running Fitch reconstruction over {n_workers} workers...")
        # One pool for the whole file (not one per block) so the spawn cost of a
        # worker is paid once; the flat tree is shipped via the initializer.
        with get_mp_context().Pool(n_workers, initializer=_init_worker, initargs=(flat,)) as pool:
            for result in pool.imap(_fitch_block_worker, blocks):
                absorb(result)

    return mutation_dict, ambiguous_total
