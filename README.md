
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
bwa
datasets --version
```

The provided `environment.yml` installs all required dependencies, including:

* Python 3.10
* BWA (classic)
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

**Note:** Multi-species mode is experimental and intended for exploratory analyses.

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

