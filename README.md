
# CORAL

**Comparative Orthologous Read-based Analysis of Lineage Substitutions**

CORAL is a tool for scalable extraction, detection, and analysis of point mutations across species evolutionary history.
It fragments each genome into overlapping pseudo-reads, aligns them to a reference genome, filters the alignments by mapping quality, extracts unambiguous trinucleotide substitutions, and summarizes mutation rates and mutation spectra.
In pair mode, the differences between two genomes are folded into 52 unpolarized trinucleotide classes. In trio mode, an outgroup is used to assign each substitution by parsimony to one of two sister lineages, giving the standard 96-category spectrum of each branch.

---

## Pipeline overview

<img width="4066" height="1176" alt="CORAL_pipeline" src="https://github.com/user-attachments/assets/dd9d9d43-8775-4585-9be7-1f0bafebfc92" />

---

## Installation

### Requirements

* Linux (or WSL2 for windows)
* Conda (Miniforge or Anaconda recommended)

### Recommended installation

```bash
git clone https://github.com/MaruvkaLab/CORAL.git
cd CORAL
conda env create -f environment.yml
conda activate coral-env
pip install -e .
```

`environment.yml` installs Python 3.10, BWA (the default aligner, BWA-MEM), SAMtools, the NCBI Datasets CLI and the
Python dependencies; `coral --help` checks the installation. Other aligners (BWA-MEM2, minimap2, BBMap) work if
installed.

### Optional: PHYLIP and repeat masking

PHYLIP is needed only for `coral run_multi` and `coral run_phylip` (`conda install -c bioconda phylip`).

Repeat masking is **off by default**. Enable it with `--repeat-mask`, which skips pseudo-reads
that are mostly repeat and suppresses calls at masked reference positions:

```bash
conda install -c bioconda blast          # windowmasker, the default backend
conda install -c bioconda repeatmasker   # only for --repeat-mask repeatmasker
```

WindowMasker needs no repeat library, so it works for any species; RepeatMasker is more precise
where a curated library exists for the clade (`--repeat-species`).

---

## Quick start

### Three-taxon pipeline (outgroup + two ingroups)

```bash
coral run_trio \
  --outgroup Saccharomyces_mikatae_IFO_1815 GCF_947241705.1 \
  --species Saccharomyces_paradoxus GCF_002079055.1 \
            Saccharomyces_cerevisiae_S288C GCF_000146045.2 \
  --output ../test_output \
  --mapq 60 \
  --suffix test
```

This runs the full pipeline, including genome download, reference indexing, pseudo-read fragmentation, alignment, mutation extraction, and summary table and plot generation.

Names should contain no spaces; double underscores separate the names in output identifiers.

---

### Two-taxon pipeline (reference + target)

```bash
coral run_pair \
  --reference Saccharomyces_mikatae_IFO_1815 GCF_947241705.1 \
  --target Saccharomyces_cerevisiae_S288C GCF_000146045.2 \
  --output ../test_output \
  --mapq 60
```

With no outgroup there is nothing to tell which lineage changed, so a call is a difference
between the two genomes rather than a mutation on a branch. Calls are folded to the 52
undirected classes (`A[C-T]G`), with the change as seen from the reference kept beside them.
The alignment and window options are those of `run_trio`; `--five-mer` is not available. The normalized
table counts each callable window half to each genome's trinucleotide.

---

### Multi-species analysis (experimental)

```bash
coral run_multi \
  --species-list '[["Drosophila_melanogaster","GCF_000001215.4"],["Drosophila_sechellia","GCF_004382195.2"],["Drosophila_mauritiana","GCF_004382145.1"],["Drosophila_simulans","GCF_016746395.2"]]' \
  --outgroup Drosophila_simulans \
  --output ../test_output \
  --run-id drosophila_test \
  --mapq 60
```

`run_multi` runs its alignment, scan and Fitch stages in parallel by default,
using every core the process was granted. To control this:

```bash
--cores N          # total core budget (default: all available)
--no-parallel      # run every stage serially
--align-jobs N     # species aligned concurrently
--scan-jobs N      # chromosomes scanned concurrently; 1 uses the whole-genome pileup
--fitch-jobs N     # workers for the Fitch pass
--max-memory-mb N  # memory budget; caps the worker count of each stage
```

`--align-jobs`, `--scan-jobs` and `--fitch-jobs` divide the `--cores` budget and cannot exceed it.

**Note:** Multi-species mode is experimental and intended for exploratory analyses.

**Note:** Parallel stages use the `spawn` start method, which re-imports
`__main__`. Driving the pipeline from the `coral` command or from a script is
fine; calling it from `python -c` or from stdin is not.

---

## Functional workflow

### Step 1: Genome preparation

* Download genomes by NCBI assembly accession
* Index the reference genome for alignment

### Step 2: Pseudo-reads and alignment

* Fragment the query genomes into 150 bp pseudo-reads with a 75 bp step (~2x tiling coverage)
* Align the pseudo-reads to the reference (the outgroup in trio mode) with BWA-MEM
* Exclude secondary alignments and alignments to unplaced or alternative contigs
* Keep pseudo-reads with MAPQ 60 (`--mapq`), and pseudo-reads with 1 ≤ MAPQ < 60 (`--low-mapq`) only if they overlap
  a neighbouring pseudo-read of the same genome on the reference (the continuity rule); MAPQ 0 pseudo-reads are excluded
* Allow customization of aligner and parameters

**Aligner constraint.** A custom aligner (`--aligner-cmd`) must write one primary record per read, in the same order as the input FASTQ.
The continuity filter treats consecutive records as neighbouring fragments of the source genome, so an aligner that reorders its output silently breaks the filter.
BWA, BWA-MEM2 and minimap2 preserve input order by default.
BBMap does not, so its command must include `ordered=t`.
The reference is indexed once per run; a custom aligner gives its index command with `--index-cmd` (placeholders `{ref}`, `{refdir}` and `{cores}`).

**MAPQ scale.** The default `--mapq 60` follows the BWA MAPQ scale. Other aligners may use different scales (e.g., BBMap never reports 60), so adjust `--mapq` and `--low-mapq` accordingly.

**Reruns.** Finished stages are reused if their output files exist, whatever the options. When changing options, use a new output directory or `--no-cache`.

### Step 3: Mutation detection

* Parse the filtered alignments along the reference
* A site is callable when it is covered by at least one retained pseudo-read, all reads agree on the base and none
  carries an insertion or deletion there (`--indel-window K` also excludes sites within K bp of an indel)
* Call substitutions only where the flanking bases are identical in all genomes
* Trio mode assigns a substitution to the sister whose base differs from the shared base of the other sister and the
  outgroup; sites where all three genomes differ are not assigned
* Optionally retain genomic positions

### Step 4: Normalization and analysis

* Normalize mutation counts by the trinucleotide abundance at callable sites (counted on the outgroup in trio mode)
* Collapse complementary strands into canonical spectra
* Generate summary tables and visualizations

---

## Output overview

Each run produces a self-contained output directory containing:

* `Mutations/*_mutations.csv.gz` - per-branch mutation lists
* `Mutations/*_mutations.json` - trinucleotide mutation counts
* `Tables/*.tsv` - normalized mutation spectra
* `Plots/*.png` - diagnostic and summary plots
* `run_summary.json` - extraction diagnostics: covered positions and why they were rejected, windows skipped over coverage gaps, and the classification of every scored site

Mutation files are named `<taxon1>__<taxon2>__<reference>__mutations.*`: the mutations inferred on the branch leading
to `taxon1` since its divergence from `taxon2`, using `reference` as the outgroup genome. Sites where both sisters
share a base that differs from the reference (`ref_differs`) cannot be polarized with three taxa, so they are not
assigned to a branch. See `OUTPUT_FORMAT.md` for the full file formats.

---

## Documentation

* `tutorial.ipynb` - command-line tutorial and examples
* `tutorial_step_by_step.ipynb` - Python API example
* `OUTPUT_FORMAT.md` - output file structure and naming conventions

---

## Citation

Details, benchmarking, and results are available in the preprint: https://doi.org/10.64898/2026.02.02.703326

The final reference will be updated upon publication.

---

