
# CORAL

**Comparative Orthologous Read-based Analysis of Lineage Substitutions**

CORAL is a tool for scalable extraction, detection, and analysis of point mutations across species evolutionary history.
It aligns multiple species to a shared reference genome, simulates reads, filters alignments by mapping quality, extracts unambiguous trinucleotide substitutions, and summarizes mutation rates and mutation spectra.

---

## Reference

Preprint available at https://doi.org/10.64898/2026.02.02.703326

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

### Verify installation

```bash
coral --help
samtools --version
bwa-mem2 version
datasets --version
```

The provided `environment.yml` installs all required dependencies, including:

* Python 3.10
* BWA-MEM2 (default aligner; `--aligner-name bwa` still works if classic BWA is installed)
* SAMtools
* NCBI Datasets CLI
* unzip
* All required Python dependencies

### Optional: PHYLIP (for phylogenetic inference)

PHYLIP is **not required** for the core pipeline.

Install only if using phylogenetic inference via `coral run_multi` or `coral run_phylip`:

```bash
conda install -c bioconda phylip
```

### Optional: repeat masking

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
coral run_single \
  --outgroup Saccharomyces_mikatae_IFO_1815 GCF_947241705.1 \
  --species Saccharomyces_paradoxus GCF_002079055.1 \
            Saccharomyces_cerevisiae_S288C GCF_000146045.2 \
  --output ../test_output \
  --mapq 60 \
  --suffix test
```

This runs the full pipeline, including genome download, reference indexing, read simulation, alignment, mutation extraction, and summary table and plot generation.

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

Setting a stage to `1`, or passing `--cores 1`, restores that stage's original
serial path. `--align-jobs`, `--scan-jobs` and `--fitch-jobs` divide the
`--cores` budget and cannot exceed it.

**Note:** Multi-species mode is experimental and intended for exploratory analyses.

**Note:** Parallel stages use the `spawn` start method, which re-imports
`__main__`. Driving the pipeline from the `coral` command or from a script is
fine; calling it from `python -c` or from stdin is not.

---

## Functional workflow

### Step 1: Genome preparation

* Download genomes by NCBI assembly accession
* Index the reference genome for alignment

### Step 2: Read simulation and alignment

* Simulate FASTQ reads by sliding a window across genomes
* Align simulated reads to the outgroup reference
* Filter alignments by MAPQ and coverage
* Allow customization of aligner and parameters

### Step 3: Mutation detection

* Generate pileups from reference and aligned BAMs
* Extract unambiguous trinucleotide substitutions
* Optionally retain genomic positions

### Step 4: Normalization and analysis

* Normalize mutation counts by underlying trinucleotide abundance
* Collapse complementary strands into canonical spectra
* Generate summary tables and visualizations

---

## Output overview

Each run produces a self-contained output directory containing:

* `Mutations/*_mutations.csv.gz` – per-branch mutation lists
* `Mutations/*_mutations.json` – trinucleotide mutation counts
* `Tables/*.tsv` – normalized mutation spectra
* `Plots/*.png` – diagnostic and summary plots
* `run_summary.json` – extraction diagnostics: pileup lines and why they were rejected, windows skipped over coverage gaps, and the classification of every scored site

`run_summary.json` also records `ref_differs` – sites where both sister taxa
share a base that differs from the reference. This is **not** a reference-branch mutation:
with only three taxa, a change on the reference branch and a change on the branch ancestral
to both sisters are equally parsimonious, so the site cannot be polarized. The accompanying
`reference_difference_spectrum` is written as `sister_base > reference_base` by convention
and mixes both directions; it is a divergence spectrum, not a branch mutation spectrum. See
`OUTPUT_FORMAT.md`.

Mutation files are named:

```
<taxon1>__<taxon2>__<reference>__mutations.*
```

This indicates mutations inferred on the branch leading to `taxon1` since divergence from `taxon2`, using `reference` as the outgroup genome.

See `OUTPUT_FORMAT.md` for full file format and naming conventions.

---

## Documentation

* `tutorial.ipynb` – command-line tutorial and examples
* `OUTPUT_FORMAT.md` – output file structure and naming conventions

---

## Citation

Details, benchmarking, and results are available in the preprint: https://doi.org/10.64898/2026.02.02.703326

The final reference will be updated upon publication.

---

