# CORAL Output Format and Naming Conventions

This document describes the output file structure and naming conventions used by CORAL.

## Output Directory Structure

```
<output_dir>/
  └── <run_id>/
      ├── Mutations/                 # Mutation lists (CSV) and counts (JSON)
      ├── Triplets/                  # Trinucleotide context files
      ├── Tables/                    # Normalized spectra tables
      ├── Plots/                     # Visualization plots
      ├── Intervals/                 # Read interval files (for coverage plots)
      ├── run_summary.json           # Extraction diagnostics
      └── pipeline_timings.json      # Pipeline execution timing information
```

The genome folders (FASTA and aligner index), `BAMs/` and any pileup are removed at the end of a run; `--no-cleanup`
keeps them. Optional outputs: `Annotated/` (`--annotate`), `Five_mers/` (`--five-mer`) and `near_indel.json`
(windows left out near indels, with `--indel-window` K > 1). If files already exist and `--no-cache` is not
specified, CORAL skips regeneration and uses the existing files.

## Run ID

The `run_id` names the output directory:

- `coral run_trio`: `<outgroup_name>__<species1_name>__<species2_name>`, plus `_<suffix>` with `--suffix`
- `coral run_pair`: `<reference_name>__<target_name>`
- `coral run_multi`: `multi_species_run`, or the value of `--run-id`

Species names use underscores (`_`) instead of spaces; double underscores (`__`) separate species names.

## Mutation Files

**Each file contains the mutations inferred to have occurred on the phylogenetic branch leading to the first listed
taxon (`<taxon1>`) since its divergence from `<taxon2>`, using `<reference>` as the reference genome.**

```
Mutations/<taxon1>__<taxon2>__<reference>__mutations.csv.gz
Mutations/<taxon1>__<taxon2>__<reference>__mutations.json
```

A trio writes one pair of files per branch. For example,
`Saccharomyces_paradoxus__Saccharomyces_cerevisiae_S288C__Saccharomyces_mikatae_IFO_1815__mutations.csv.gz` contains
the mutations on the branch leading to `Saccharomyces_paradoxus` after its split from `Saccharomyces_cerevisiae_S288C`,
with `Saccharomyces_mikatae_IFO_1815` as the outgroup.

The CSV has one row per mutation: `chromosome`, `position` (1-based, on the reference) and `mutation` (e.g.
`A[C>T]G`). The JSON counts the mutations by trinucleotide substitution type, before strand collapse.

### Pair mode

`coral run_pair` writes one pair of files, for the target against the reference:
```
<target>__<reference>__mutations.csv.gz
<target>__<reference>__mutations.json
```

The CSV has a row per call: `chromosome,position,change,mutation`, for example
`24,10976,C[T>C]T,C[C-T]T`. `change` is what the window looks like from the reference
(reference `CTT`, target `CCT`); `mutation` is the same call folded to one of the 52
undirected classes. The JSON counts those 52 classes. Because there is no outgroup, neither
column says which lineage changed.

## Triplet Files

`Triplets/<taxon1>__<taxon2>__<reference>__triplets.json` counts the outgroup trinucleotides at the branch's callable
sites (identical in all three species or mutated in one sister), the denominator of the normalization.

Pair mode writes `<target>__<reference>__triplets.json`: every window that was called or
identical, i.e. the windows a call could have been made in, counted half to the reference 3-mer
and half to the target 3-mer (identical windows count 1, a called window 0.5 to each allele).

## Tables

Trio mode writes, with one column per mutation file:

- `normalized_scaled.tsv` - Normalized mutation spectra (scaled)
- `collapsed_mutations.tsv` - Raw mutation counts (collapsed)
- `scaled_raw.tsv` - Scaled raw mutation counts
- `triplets.tsv` - Trinucleotide context counts

Rows are the 96 mutation types `X[Y>Z]W` (`X` and `W` the flanking bases, `Y` the reference base, `Z` the derived
base; `A[T>C]A` is `ATA` to `ACA`), or the 32 trinucleotide contexts in `triplets.tsv`. For example:

```
	Saccharomyces_cerevisiae_S288C__..._mutations	Saccharomyces_paradoxus__..._mutations
A[T>C]A	4812	3749
T[T>C]T	6118	5466
T[T>A]C	352	236
```

Pair mode writes `Tables/<target>__<reference>__normalized.tsv`, with the columns `label`, `count`,
`opportunity`, `rate` (count / opportunity) and `normalized` (rates scaled to sum to 10,000). A class's opportunity
counts both of its trinucleotides, on either strand.

## Intervals and Plots

`Intervals/<species>_to_<reference>_intervals.tsv.gz` lists the reference intervals covered by aligned reads
(`chromosome`, `start`, `end`; 0-based, end-exclusive), used for the coverage plots.

`Plots/` holds the mutation spectra (`*_normalized.png`, `*_raw.png`, `*_triplets.png`), coverage
(`coverage_<chromosome>.png`) and mutation density (`mutation_density_<chromosome>.png`) plots. The MAPQ histograms
are written to `<output_dir>/Plots/<species>_to_<reference>.png`.

## Multi-Species Pipeline Files

For `coral run_multi`, additional files are created in `<run_id>/`:

- `matching_bases.csv.gz` - Mutation matrix for phylogenetic analysis
- `annotated_tree.nwk` - Newick tree with branch annotations
- `species_mapping.json` - Mapping between species names and internal IDs
- `triplets.json` - Trinucleotide context counts from the scan, used as the
  normalization denominator
- `mutation_spectras.tsv` - Mutation spectra summary

`run_multi` also writes `Mutations/` and `Tables/` in the formats described
above, one column per tree branch instead of one per species pair. Branch names
come from PHYLIP's interior numbering, so they match the `.outfile` of the
corresponding run in `multi_species_phylip_no_tree/` or
`multi_species_phylip_with_tree/`.

## Run summary (`run_summary.json`)

One file per run, written by the mutation extractor. It records
what the scan saw, including the sites that produce no mutation call and therefore appear
in no other output. Counts are omitted when zero, so read them with a default of 0.

```json
{
  "pileup_lines":      { "lines_after_mask": 4166311, "lines_masked": 0,
                         "no_depth": 12034, "reads_disagree": 881, "insertion": 55 },
  "candidate_windows": { "non_consecutive": 20714 },
  "site_classes":      { "identical": 2564353, "taxa1_mut": 105862, "taxa2_mut": 136417,
                         "ref_differs": 227255, "all_differ": 22913,
                         "flanks_not_conserved": 1109511 },
  "reference_difference_spectrum": { "A[C>T]G": 3311, "...": 0 }
}
```

**`pileup_lines`** counts covered reference positions (one pileup line each). `lines_after_mask` is the number of lines
the scan actually saw, i.e. after repeat masking removed any; `lines_masked` is how many
masking removed. The remaining keys are per-line QC rejections: `no_depth`, `deletion`,
`insertion`, `reads_disagree`, `unparsed`. A line is attributed to the *first* check it
fails, so a line with both an insertion and low depth is counted once, under `insertion`.
`reads_disagree` means reads within one species disagree at that position; pseudo-reads
tiling one assembly are identical where they overlap, so this indicates reads arriving from
different loci (paralogy or mismapping).

**`candidate_windows`** counts 3-position windows, not lines. `non_consecutive` is windows
whose three lines pass QC and lie on one contig but are not at adjacent positions, i.e. the
window spans a coverage gap. Windows spanning a contig boundary are never formed and are
not counted, so that the whole-genome scan and the per-chromosome parallel scan produce
identical summaries.
`near_indel` (only with `--indel-window` K > 1) is windows left out because a read has a
deleted base or an insertion within K bp of their middle base.

**`site_classes`** in pair mode are `differs` (a call), `identical`, `flanks_not_conserved`
and `not_acgt` (a window base is not A/C/G/T in the reference or the target).

**`site_classes`** counts the windows that were scored, partitioning them:

| Reference | Taxa1 | Taxa2 | Class |
|---|---|---|---|
| A | A | A | `identical` |
| A | C | A | `taxa1_mut` |
| A | A | C | `taxa2_mut` |
| A | C | C | `ref_differs` |
| A | C | G | `all_differ` |

with `flanks_not_conserved` for any window whose flanking bases are not identical across
all three species. `identical + taxa1_mut + taxa2_mut` is the callable-site total that the
triplet counts and mutation rates are based on.

**`ref_differs` is deliberately not called an outgroup mutation.** Both
sisters carry one base and the reference carries another, but three taxa cannot polarize
that: a single change on the reference branch and a single change on the branch ancestral
to both sisters are equally parsimonious. Resolving it needs a fourth taxon.

**`reference_difference_spectrum`** is the trinucleotide spectrum of those same sites,
written as `sister_base > reference_base`. Because the change is unpolarized, this
direction is a convention, not an inference - where the change actually occurred on the
sister-ancestor branch, the true substitution is the reverse. Treat it as a divergence
spectrum between the sister pair's ancestor and the reference, not as a branch-specific
mutation spectrum, and do not feed it to signature extraction as if it were one.

The pipeline adds four more sections after extraction: **`run`** (run id, taxa, parameters,
tool versions), **`genomes`** (accession, contig count and total bp), **`alignment`**
(per species: read counts split into `kept_high_mapq` versus `kept_rescued` by the
continuity rule, the filtered categories, and the full MAPQ histogram), and **`timings`**
(per-stage wall and CPU seconds, per-species alignment time).
`alignment` is absent when the BAMs were reused from a previous run rather than re-filtered.
