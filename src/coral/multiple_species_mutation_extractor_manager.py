import csv
import os
import gzip
from collections import defaultdict
import pandas as pd
import json
from .multiple_species_utils import (annotate_tree_with_indices, 
                                     save_annotated_tree, 
                                     collapse_mutations, 
                                     filter_mutations_dict,
                                     collapse_triplets,
                                     filter_triplets_dict,
                                     normalize_by_triplets,
                                     scale_counts)
from .plot_utils import MutationSpectraPlotter
from .utils import log
import re


MIN_DEPTH = 1  # Minimum depth threshold for quality check

# Matches read-start (^ + mapQ char), read-end ($), or an indel marker (+N / -N).
# Digits are captured so the N following indel bases can be skipped in one pass.
_CLEAN_RE = re.compile(r'\^.|\$|[+-](\d*)')


def clean_bases(s):
    """Strip read-start/end markers and indel notation from an mpileup bases field.

    Fast path: every marker this removes (^x, $, +N…, -N…) contains one of the
    four characters ^ $ + - . When none are present — the overwhelmingly common
    case for CORAL's shallow pileups — the regex below matches nothing and would
    return s unchanged, so we skip it entirely. Byte-for-byte identical to the
    regex path but avoids the regex-engine call and list building.
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
    """Detects the multi-species sites and runs Fitch over them.

    Two stages have both a serial and a parallel implementation, and which runs
    is decided per stage:

    * **the pileup scan** -- serially, one Python loop over the whole-genome
      pileup at ``pileup_file``; in parallel, one ``samtools mpileup -r`` per
      chromosome straight from ``bams``, so ``pileup_file`` is never needed.
      Enabled by passing ``ref_fasta``, ``bams`` and ``fai_path``.
    * **the Fitch pass** -- serially, the recursive walk with an ete3
      ``tree.copy()`` per site; in parallel, the flattened, copy-free walk over
      row blocks. Enabled with ``parallel_fitch=True``.

    Both parallel paths are meant to produce output identical to the serial ones;
    ``tests/`` checks that by running the two against each other. Constructed
    with none of the parallel arguments, this class behaves exactly as it always
    has.
    """

    def __init__(self, pileup_file, output_dir, n_species, tree=None, species_list=None, mapping=None, no_cache=False, verbose=False,
                 ref_fasta=None, bams=None, fai_path=None, cores=1, scan_jobs=None, fitch_jobs=None, max_memory_mb=None, parallel_fitch=False):
        self.pileup_file = pileup_file
        self.output_dir = output_dir
        self.n_species = n_species
        self.tree = tree
        self.species_list = species_list
        self.mapping = mapping
        self.no_cache = no_cache
        self.verbose = verbose

        # -- parallel configuration (all optional; defaults reproduce the
        # -- original serial behaviour exactly)
        self.ref_fasta = ref_fasta
        self.bams = list(bams) if bams else None
        self.fai_path = fai_path
        self.cores = cores or 1
        self.scan_jobs = scan_jobs
        self.fitch_jobs = fitch_jobs
        self.max_memory_mb = max_memory_mb
        self.parallel_fitch = parallel_fitch
        if self.tree is None and self.species_list is None:
            raise ValueError("Either newick_tree or species_list must be provided.")
        if self.mapping is None:
            raise ValueError("Dictionary mapping taxa names must be provided.")

        os.makedirs(self.output_dir, exist_ok=True)
        self.plots_dir = os.path.join(self.output_dir, "Plots")
        self.csv_dir = os.path.join(self.output_dir, "Mutations")

        self.triplet_counts = {}   # genome-wide, branch-agnostic triplet denominator

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

        A flank-clean site (every sample equals its own line's ref base at prev
        and next) falls into one of three center cases, named to match the
        two-taxa detect_mutation_triplet():

          (1) ingroup varies among themselves -> emit a CSV row for Fitch AND
              count the context (the multi-allelic-center case you accepted).
          (2) every taxon INCLUDING the outgroup/ref (taxa0) shows the ref base
              -> no row; count the context (a zero-mutation opportunity).
          (3) the ingroup is uniform but differs from the outgroup/ref base
              -> a substitution on the ingroup-ancestor branch that a single
              outgroup cannot polarise. This is exactly the two-taxa
              "both taxa differ from ref" case, which the two-taxa triplet
              counter EXCLUDES -> not a usable opportunity: no row and NOT
              counted in the denominator. (This fixes limitation 2.)

        matches_ref/normalize are computed once and shared by all three cases.
        """
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

        if not (matches_ref(buffer[0]) and matches_ref(buffer[2])):
            return None, None

        ref_upper = buffer[1][2].upper()
        curr_bases = normalize(buffer[1], buffer[1][2])   # already uppercased
        distinct = set(curr_bases)

        # Case (3): ingroup uniform but != outgroup/ref -> excluded everywhere.
        # count a triplet opportunity only when >=1 species matches the ref at the
        # center -> drops BOTH a==b!=r (case 3) and the tri-allelic a!=b,both!=r case,
        # exactly matching two-taxa detect_mutation_triplet.
        if ref_upper not in distinct: # if len(distinct) == 1 and ref_upper not in distinct:
            return None, None
        triplet_context = buffer[0][2].upper() + ref_upper + buffer[2][2].upper()


        if len(distinct) > 1:   # case (1): variable across ingroup -> CSV row
            mutation_row = [
                buffer[1][0],           # chrom
                buffer[1][1],           # pos
                buffer[0][2].upper(),   # left  = prev line's ref base
                buffer[2][2].upper(),   # right = next line's ref base
                ref_upper,
            ] + [b.upper() for b in curr_bases]
            return mutation_row, triplet_context

        return None, triplet_context   # case (2): fully invariant


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

    @property
    def parallel_scan(self):
        """True when the scan can run per chromosome: it needs the indexed BAMs
        and the reference, which the serial path does not."""
        return bool(self.ref_fasta and self.bams and self.fai_path)

    def _scan_pileup_serial(self, csv_path, triplets_path, header):
        """The original scan: one Python loop over the whole-genome pileup.

        Kept verbatim. It is both the default path and the oracle the parallel
        scan is tested against, so it should not be "tidied" alongside it.
        """
        triplet_counts = defaultdict(int)

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
                        result, triplet = self._detect_site(buffer)
                        if triplet is not None:
                            triplet_counts[triplet] += 1
                        if result is not None:
                            writer.writerow(result)

        triplet_counts = dict(triplet_counts)
        with open(triplets_path, 'w') as tf:
            json.dump(triplet_counts, tf, indent=2)
        log(f"Saved {len(triplet_counts)} triplet contexts to {triplets_path}", self.verbose)
        return triplet_counts

    def _fitch_serial(self, csv_path):
        """The original Fitch pass: recursive, with an ete3 deep copy per site."""
        mutation_dict = defaultdict(list)
        ambiguous_counter = 0

        for chunk in pd.read_csv(csv_path, chunksize=1000):
            for _, row in chunk.iterrows():
                mutation_dict, ambiguous = self._fitch(self.tree.copy(), row, mutation_dict)
                ambiguous_counter += ambiguous

        return mutation_dict, ambiguous_counter

    def extract(self):
        csv_path = os.path.join(self.output_dir, "matching_bases.csv.gz")
        triplets_path = os.path.join(self.output_dir, "triplets.json")


        header = ["chromosome", "position", "left", "right"] + [f"taxa{k}" for k in self.mapping if isinstance(k, int)]

        
        cache_ok = (not self.no_cache
                    and os.path.exists(csv_path)
                    and os.path.exists(triplets_path))
        
        if cache_ok:
            log(f'Using cached matching positions from csv at {csv_path}', self.verbose)
            with open(triplets_path) as tf:
                self.triplet_counts = json.load(tf)
        elif self.parallel_scan:
            # Imported here so the serial path never pays for multiprocessing.
            from .parallel_scan import scan_pileup_parallel
            self.triplet_counts = scan_pileup_parallel(self, csv_path, triplets_path, header)
        else:
            self.triplet_counts = self._scan_pileup_serial(csv_path, triplets_path, header)

        if self.tree:
            if self.parallel_fitch:
                from .parallel_fitch import run_fitch
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

    def _save_results(self, mutation_dict):
        spectra_plotter = MutationSpectraPlotter()
        os.makedirs(self.plots_dir, exist_ok=True)
        os.makedirs(self.csv_dir, exist_ok=True)
        tables_dir = os.path.join(self.output_dir, "Tables")
        os.makedirs(tables_dir, exist_ok=True)

        # Single shared, branch-agnostic denominator.
        collapsed_triplets = collapse_triplets(filter_triplets_dict(self.triplet_counts or {}))

        spectra_dict = {}       # per-branch collapsed+filtered raw spectrum
        normalized_scaled = {}  # per-branch (spectrum / shared triplet vector), scaled
        scaled_raw = {}         # per-branch raw spectrum, scaled

        for branch_key, mutations in mutation_dict.items():
            df = pd.DataFrame(mutations, columns=["chromosome", "position", "mutation"])
            csv_path = os.path.join(self.csv_dir, f"{branch_key}.csv.gz")
            df.to_csv(csv_path, index=False, header=False, sep="\t", compression="gzip")
            mutation_spectra = collapse_mutations(dict(df['mutation'].value_counts()))
            mutation_spectra = filter_mutations_dict(mutation_spectra)
            spectra_dict[branch_key] = mutation_spectra

            normalized_scaled[branch_key] = scale_counts(normalize_by_triplets(mutation_spectra, collapsed_triplets))
            scaled_raw[branch_key] = scale_counts(mutation_spectra)
            
            spectra_plot_path = os.path.join(self.plots_dir, f"{branch_key}_spectra.png")
            spectra_plotter.plot_mutations(pd.Series(mutation_spectra), spectra_plot_path, f"Mutation Spectra: {branch_key}")
            
        spectra_df = pd.DataFrame(spectra_dict)
        spectra_df.to_csv(os.path.join(self.output_dir, "mutation_spectras.tsv"), sep="\t")

        pd.DataFrame(spectra_dict).to_csv(os.path.join(tables_dir, "collapsed_mutations.tsv"), sep="\t")
        pd.DataFrame(normalized_scaled).to_csv(os.path.join(tables_dir, "normalized_scaled.tsv"), sep="\t")
        pd.DataFrame(scaled_raw).to_csv(os.path.join(tables_dir, "scaled_raw.tsv"), sep="\t")
        # Denominator is one shared vector -> single column (not per-branch).
        pd.Series(collapsed_triplets, name="triplets").to_frame().to_csv(os.path.join(tables_dir, "triplets.tsv"), sep="\t")

        
