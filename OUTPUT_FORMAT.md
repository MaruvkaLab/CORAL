# CORAL Output Format and Naming Conventions

This document describes the output file structure and naming conventions used by CORAL.

## Output Directory Structure

CORAL creates a self-contained output directory for each run. The directory structure follows this pattern:

```
<output_dir>/
  └── <run_id>/
      ├── *.fasta                    # Genome FASTA files
      ├── *.pileup.gz               # Multi-taxa pileup file
      ├── *_mutations.csv.gz        # Full mutation lists (one per species pair)
      ├── *_mutations.json          # Mutation context counts (one per species pair)
      ├── Mutations/                 # Mutation files directory
      ├── Triplets/                  # Trinucleotide context files
      ├── Tables/                    # Normalized spectra tables
      ├── Plots/                     # Visualization plots
      ├── Intervals/                 # Read interval files (for coverage plots)
      ├── run_summary.json           # Extraction diagnostics (one per run/trio)
      └── pipeline_timings.json      # Pipeline execution timing information
```

## Run ID Naming Convention

The `run_id` determines the output directory name and is used in many file names.

### For `coral run_single`:

**Default run_id:**
```
<outgroup_name>__<species1_name>__<species2_name>
```

**Example:**
- Outgroup: `Saccharomyces_mikatae_IFO_1815`
- Species 1: `Saccharomyces_paradoxus`
- Species 2: `Saccharomyces_cerevisiae_S288C`
- **Run ID:** `Saccharomyces_mikatae_IFO_1815__Saccharomyces_paradoxus__Saccharomyces_cerevisiae_S288C`

**With `--suffix`:**
If you specify `--suffix test`, the run_id becomes:
```
<outgroup_name>__<species1_name>__<species2_name>_test
```

### For `coral run_multi`:

**Default run_id:**
```
multi_species_run
```

**Custom run_id:**
If you specify `--run-id <custom_id>`, that value is used directly.

## File Naming Conventions

### Pileup Files

**Pattern:**
```
<run_id>.pileup.gz
```

**Example:**
```
Saccharomyces_mikatae_IFO_1815__Saccharomyces_paradoxus__Saccharomyces_cerevisiae_S288C.pileup.gz
```

### Mutation Files

Mutation files are named based on a species pair and a reference genome. **Each file contains the mutations inferred to have occurred on the phylogenetic branch leading to the first listed taxon (`<taxon1>`) since its divergence from `<taxon2>`, using `<reference>` as the reference genome.** In other words, for the file `<taxon1>__<taxon2>__<reference>__mutations.csv.gz`, the mutations listed are those on the branch leading to `<taxon1>`, relative to the common ancestor with `<taxon2>`.

**Pattern:**
```
<taxon1>__<taxon2>__<reference>__mutations.csv.gz
<taxon1>__<taxon2>__<reference>__mutations.json
```

**Example:**
For a run with:
- Reference: `Saccharomyces_mikatae_IFO_1815`
- Taxon 1: `Saccharomyces_paradoxus`
- Taxon 2: `Saccharomyces_cerevisiae_S288C`

Two mutation files are created for each direction (i.e., branch), corresponding to the mutations on the branch leading to each taxon:
```
Saccharomyces_paradoxus__Saccharomyces_cerevisiae_S288C__Saccharomyces_mikatae_IFO_1815__mutations.csv.gz      # Mutations on the S. paradoxus branch
Saccharomyces_paradoxus__Saccharomyces_cerevisiae_S288C__Saccharomyces_mikatae_IFO_1815__mutations.json
Saccharomyces_cerevisiae_S288C__Saccharomyces_paradoxus__Saccharomyces_mikatae_IFO_1815__mutations.csv.gz      # Mutations on the S. cerevisiae branch
Saccharomyces_cerevisiae_S288C__Saccharomyces_paradoxus__Saccharomyces_mikatae_IFO_1815__mutations.json
```
*For example, `Saccharomyces_paradoxus__Saccharomyces_cerevisiae_S288C__Saccharomyces_mikatae_IFO_1815__mutations.csv.gz` contains all mutations inferred to have occurred on the branch leading to `Saccharomyces_paradoxus` after its split from `Saccharomyces_cerevisiae_S288C`, using `Saccharomyces_mikatae_IFO_1815` as the outgroup/reference.*

**Location:**
- CSV files: `<run_id>/Mutations/`
- JSON files: `<run_id>/Mutations/`

### Triplet Files

**Pattern:**
```
<taxon1>__<taxon2>__<reference>__triplets.json
```

**Location:**
`<run_id>/Triplets/`

### Normalized Spectra Tables

**Pattern:**
The following files are created in the `Tables/` directory:

- `normalized_scaled.tsv` - Normalized mutation spectra (scaled)
- `collapsed_mutations.tsv` - Raw mutation counts (collapsed)
- `scaled_raw.tsv` - Scaled raw mutation counts
- `triplets.tsv` - Trinucleotide context counts

**Location:**
`<run_id>/Tables/`

**Table Format:**

#### `triplets.tsv`

This file contains trinucleotide context counts for each mutation file. Each row represents a trinucleotide context (e.g., `ATG`, `ATT`, `ATA`), and columns represent different mutation files.

**Format:**
- First column: Trinucleotide context (e.g., `ATG`, `ATT`, `ATA`)
- Subsequent columns: Counts for each mutation file, with column headers matching the mutation file names (without extensions)

**Example:**
```
Saccharomyces_cerevisiae_S288C__Saccharomyces_paradoxus__Saccharomyces_mikatae_IFO_1815__mutations	Saccharomyces_paradoxus__Saccharomyces_cerevisiae_S288C__Saccharomyces_mikatae_IFO_1815__mutations
ATG	120647	120647
ATT	165537	165537
ATA	128936	128936
TTC	138995	138995
CCT	56896	56896
TTT	213608	213608
```

#### `normalized_scaled.tsv`

This file contains normalized mutation spectra with mutation types in the format `X[Y>Z]W`, where:
- `X` and `W` are the 5' and 3' flanking bases (trinucleotide context)
- `Y` is the reference base
- `Z` is the derived base

Each row represents a mutation type, and columns represent counts for each mutation file.

**Format:**
- First column: Mutation type in format `X[Y>Z]W` (e.g., `A[T>C]A`, `T[T>C]T`)
- Subsequent columns: Normalized counts for each mutation file, with column headers matching the mutation file names (without extensions)

**Example:**
```
Saccharomyces_cerevisiae_S288C__Saccharomyces_paradoxus__Saccharomyces_mikatae_IFO_1815__mutations	Saccharomyces_paradoxus__Saccharomyces_cerevisiae_S288C__Saccharomyces_mikatae_IFO_1815__mutations
A[T>C]A	4812	3749
T[T>C]T	6118	5466
T[T>A]C	352	236
T[C>T]T	5322	3368
C[C>T]T	3558	2414
C[T>C]G	3635	3480
```

**Note:** The mutation type format `X[Y>Z]W` represents a substitution from `Y` to `Z` in the context of flanking bases `X` and `W`. For example, `A[T>C]A` means a T→C substitution in the context of A on both sides (i.e., the trinucleotide context is `ATA` → `ACA`).

### Interval Files

**Pattern:**
```
<base_bam_name>_intervals.tsv.gz
```

Where `<base_bam_name>` is derived from the BAM file name (without `.bam` extension).

**Example:**
```
Saccharomyces_paradoxus_to_Saccharomyces_mikatae_IFO_1815_intervals.tsv.gz
```

**Location:**
`<run_id>/Intervals/`

### Plot Files

**Pattern:**
Plots are generated with descriptive names:

- Mutation spectra: `*_normalized.png`, `*_raw.png`, `*_triplets.png`
- Coverage plots: `coverage_<chromosome>.png`
- Mutation density: `mutation_density_<chromosome>.png`
- MAPQ histograms: `<species>_to_<reference>_MAPQ.png`

**Location:**
`<run_id>/Plots/`

### Multi-Species Pipeline Files

For `coral run_multi`, additional files are created:

- `matching_bases.csv.gz` - Mutation matrix for phylogenetic analysis
- `annotated_tree.nwk` - Newick tree with branch annotations
- `species_mapping.json` - Mapping between species names and internal IDs
- `mutation_spectras.tsv` - Mutation spectra summary

**Location:**
`<run_id>/`

## File Formats

### Pileup Files (`.pileup.gz`)
- Format: Gzipped text file
- Content: Multi-taxa pileup format from samtools mpileup
- Contains: Reference and aligned species base calls at each position

### Mutation CSV Files (`.csv.gz`)
- Format: Gzipped CSV
- Columns: `chromosome`, `position`, `reference_base`, `taxon1_base`, `taxon2_base`, `context`, etc.
- Contains: Full list of detected mutations with genomic positions

### Mutation JSON Files (`.json`)
- Format: JSON
- Content: Mutation context counts (trinucleotide substitution counts)
- Structure: Dictionary mapping mutation types to counts

### TSV Tables (`.tsv`)
- Format: Tab-separated values
- Content: Normalized mutation spectra, collapsed counts, or triplet frequencies
- Structure: 
  - First column: Row identifier (trinucleotide context or mutation type)
  - Subsequent columns: Counts for each mutation file (column headers match mutation file names)
- Mutation type format: `X[Y>Z]W` where `X` and `W` are flanking bases, `Y` is reference base, `Z` is derived base
- See [Normalized Spectra Tables](#normalized-spectra-tables) section above for detailed format examples

### Interval Files (`.tsv.gz`)
- Format: Gzipped TSV
- Columns: `chromosome`, `start`, `end`
- Content: Genomic intervals covered by aligned reads

### Run Summary (`run_summary.json`)

One file per run directory, i.e. per trio, written by the mutation extractor. It records
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

**`pileup_lines`** counts individual pileup lines. `lines_after_mask` is the number of lines
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

The pipeline adds four more sections after extraction: **`run`** (run id, taxa, parameters,
tool versions), **`genomes`** (accession, contig count and total bp), **`alignment`**
(per species: read counts split into `kept_high_mapq` versus `kept_rescued` by the
continuity rule, the filtered categories, and the full MAPQ histogram), and **`timings`**
(per-stage wall and CPU seconds, per-species alignment time).
`alignment` is absent when the BAMs were reused from a previous run rather than re-filtered.

**`reference_difference_spectrum`** is the trinucleotide spectrum of those same sites,
written as `sister_base > reference_base`. Because the change is unpolarized, this
direction is a convention, not an inference — where the change actually occurred on the
sister-ancestor branch, the true substitution is the reverse. Treat it as a divergence
spectrum between the sister pair's ancestor and the reference, not as a branch-specific
mutation spectrum, and do not feed it to signature extraction as if it were one.

## Notes

1. **Species name format**: Species names use underscores (`_`) instead of spaces, matching NCBI naming conventions.

2. **Separator**: Double underscores (`__`) are used to separate species names in file names to avoid ambiguity.

3. **Compression**: Large files (pileups, mutation CSVs, intervals) are gzipped to save space.

4. **Directory organization**: Files are organized into subdirectories (`Mutations/`, `Tables/`, `Plots/`, etc.) for clarity.

5. **Caching**: If files already exist and `--no-cache` is not specified, CORAL will skip regeneration and use existing files.

