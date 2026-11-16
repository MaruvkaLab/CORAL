# Diagnostic plots: spectra, read coverage and mutation density along chromosomes

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re
from typing import List, Tuple, Optional

from .utils import log

COLOR_MUTATION = {
    "C>A": "#E64B35", "C>G": "#4DBBD5", "C>T": "#00A087",
    "T>A": "#3C5488", "T>C": "#F39B7F", "T>G": "#8491B4"
}

COLOR_TRIPLET = {
    "C": "#7E6148", "T": "#0B0B0A"
}

class MutationSpectraPlotter:
    @staticmethod
    def plot_mutations(series, output_path, title):
        df = series.reset_index()
        df.columns = ["index", "count"]
        df[['First_Base', 'Mutation', 'Third_Base']] = df['index'].str.extract(r'(\w)\[(\w>\w)\](\w)')
        df = df.sort_values(by=['Mutation', 'First_Base', 'Third_Base'])
        df.index = df["index"]
        sorted_data = df["count"]
        colors = [COLOR_MUTATION[m.split("[")[1][:-2]] for m in sorted_data.index]

        plt.figure(figsize=(12, 5), dpi=300)
        plt.bar(sorted_data.index, sorted_data.values, color=colors)
        plt.xticks(rotation=90, fontsize=6)
        plt.ylabel("Count")
        plt.title(title)
        plt.tight_layout()
        plt.savefig(output_path)
        plt.close()

    @staticmethod
    def plot_triplets(series, output_path, title):
        df = series.sort_index()
        colors = [COLOR_TRIPLET.get(t[1], "gray") for t in df.index]

        plt.figure(figsize=(12, 5), dpi=300)
        plt.bar(df.index, df.values, color=colors)
        plt.xticks(rotation=90, fontsize=6)
        plt.ylabel("Triplet Count")
        plt.title(title)
        plt.tight_layout()
        plt.savefig(output_path)
        plt.close()

    @staticmethod
    def plot_mutation_spectra_overlay(data1, data2, labels, file_name=None):
        def prepare_data(data):
            df = data.reset_index()
            df.columns = ["index", "count"]
            df[['First_Base', 'Mutation', 'Third_Base']] = df['index'].str.extract(r'(\w)\[(\w>\w)\](\w)')
            df = df.sort_values(by=['Mutation', 'First_Base', 'Third_Base'])
            df.index = df["index"]
            return df["count"]

        sorted_data1 = prepare_data(data1)
        sorted_data2 = prepare_data(data2)

        categories = sorted_data1.index
        color_dict = {"C>A": "red", "C>G": "green", "C>T": "blue", "T>A": "orange", "T>C": "purple", "T>G": "brown"}
        tick_colors = [color_dict[m.split("[")[1][:-2]] for m in categories]

        fig, ax = plt.subplots(figsize=(14, 6), dpi=300)
        x = range(len(categories))

        ax.bar(x, sorted_data2.values, color="red", width=0.6, label=labels[1], align='center', alpha=0.5)
        ax.bar(x, sorted_data1.values, color="yellow", width=0.6, label=labels[0], align='center', alpha=0.5)

        ax.set_xticks(x)
        ax.set_xticklabels(categories, fontsize=6, rotation=90)
        for tick, color in zip(ax.get_xticklabels(), tick_colors):
            tick.set_color(color)

        plt.xlabel("Mutation category")
        plt.ylabel("Mutation count")
        plt.title(f"Mutation Spectra Comparison: {labels[0]} vs {labels[1]}")
        plt.legend()
        plt.tight_layout()
        if file_name:
            plt.savefig(file_name)
        else:
            plt.show()
        plt.close()
    
    def plot(self, tables_dir, output_dir=None, verbose=True):
        """Plot raw, normalized, triplet and overlay spectra from the TSVs in tables_dir (default output: sibling Plots/)."""
        input_dir = tables_dir
        if output_dir is None:
            output_dir = os.path.join(os.path.dirname(input_dir), "Plots")
        os.makedirs(output_dir, exist_ok=True)

        def log(msg):
            if verbose:
                print(f"[plot_runner] {msg}")

        log("Loading mutation summary tables...")
        norm = pd.read_csv(os.path.join(input_dir, "normalized_scaled.tsv"), sep='\t', index_col=0)
        raw = pd.read_csv(os.path.join(input_dir, "collapsed_mutations.tsv"), sep='\t', index_col=0)
        scaled = pd.read_csv(os.path.join(input_dir, "scaled_raw.tsv"), sep='\t', index_col=0)
        trip = pd.read_csv(os.path.join(input_dir, "triplets.tsv"), sep='\t', index_col=0)

        for col in norm.columns:
            log(f"Plotting normalized mutations for {col}")
            self.plot_mutations(
                norm[col], os.path.join(output_dir, f"{col}_normalized.png"),
                f"Normalized Mutation Spectrum: {col}"
            )

        for col in raw.columns:
            log(f"Plotting raw mutations for {col}")
            self.plot_mutations(
                raw[col], os.path.join(output_dir, f"{col}_raw.png"),
                f"Raw Mutation Spectrum: {col}"
            )

        for col in trip.columns:
            log(f"Plotting triplet counts for {col}")
            self.plot_triplets(
                trip[col], os.path.join(output_dir, f"{col}_triplets.png"),
                f"Triplet Spectrum: {col}"
            )

        if len(norm.columns) == 2:
            species1, species2 = norm.columns
            log(f"Plotting overlay for normalized spectra: {species1} vs {species2}")
            self.plot_mutation_spectra_overlay(
                norm[species1], norm[species2], labels=[species1, species2],
                file_name=os.path.join(output_dir, f"{species1}_vs_{species2}_normalized_overlay.png")
            )

            log(f"Plotting overlay for scaled spectra: {species1} vs {species2}")
            self.plot_mutation_spectra_overlay(
                scaled[species1], scaled[species2], labels=[species1, species2],
                file_name=os.path.join(output_dir, f"{species1}_vs_{species2}_overlay.png")
            )



class CoveragePlotter:
    def __init__(self, fai_file, verbose = True):
        self.chrom_lengths = self._parse_fai(fai_file)
        self.verbose = verbose

    def _parse_fai(self, fai_file):
        chrom_lengths = {}
        with open(fai_file) as f:
            for line in f:
                fields = line.strip().split('\t')
                chrom = fields[0]
                length = int(fields[1])
                chrom_lengths[chrom] = length
        return chrom_lengths

    def compute_binned_coverage(self, interval_file, chrom, bin_size=1000, slide=1000):
        """Return bin midpoints and mean depth per bin (summed interval overlap / bin_size)."""
        chrom_length = self.chrom_lengths.get(chrom)
        if chrom_length is None:
            raise ValueError(f"Chromosome {chrom} not found in .fai index.")

        df = pd.read_csv(
            interval_file,
            sep='\t',
            compression='infer',
            dtype={"chromosome": str, "start": int, "end": int},
            header=0
        )
        df = df[df["chromosome"] == chrom].sort_values("start").reset_index(drop=True)

        starts = np.arange(0, chrom_length - bin_size + 1, slide)
        coverage = []
        current_index = 0
        n = len(df)

        for bin_start in starts:
            bin_end = bin_start + bin_size
            total_overlap = 0
            while current_index < n and df.at[current_index, "end"] <= bin_start:
                current_index += 1
            check_index = current_index
            while check_index < n and df.at[check_index, "start"] < bin_end:
                read_start = df.at[check_index, "start"]
                read_end = df.at[check_index, "end"]
                overlap_start = max(read_start, bin_start)
                overlap_end = min(read_end, bin_end)
                total_overlap += max(0, overlap_end - overlap_start)
                check_index += 1
            coverage.append(total_overlap / bin_size)

        midpoints = [start + bin_size // 2 for start in starts]
        return midpoints, coverage

    def compute_coverage_for_normalization(self, interval_file, chrom, bin_size=1000, slide=1000):
        """Return only the per-bin depths, for normalizing mutation counts."""
        _, coverage = self.compute_binned_coverage(interval_file, chrom, bin_size, slide)
        return coverage

    def plot_coverage(self, midpoints_list, coverage_list, labels, chrom, output_path):
        plt.figure(figsize=(15, 5))
        for midpoints, coverage, label in zip(midpoints_list, coverage_list, labels):
            plt.plot(midpoints, coverage, label=label, lw=1.5)

        plt.title(f"Coverage over {chrom} — {' vs '.join(labels)}")
        plt.xlabel("Genomic Position")
        plt.ylabel("Normalized Read Coverage (per base)")
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(output_path, dpi=300)
        plt.close()

        log(f"Saved coverage plot: {output_path}", self.verbose)

    def plot(
        self,
        interval_dir: str,
        chromosome: str,
        output_dir: str,
        bin_size: int = 100000,
        slide: Optional[int] = None
    ):
        if slide is None:
            slide = bin_size

        if chromosome not in self.chrom_lengths:
            raise ValueError(f"Chromosome {chromosome} not found in FAI file.")

        interval_files = sorted([
            os.path.join(interval_dir, f)
            for f in os.listdir(interval_dir)
            if f.endswith(".tsv") or f.endswith(".tsv.gz")
        ])
        labels = [os.path.splitext(os.path.splitext(os.path.basename(f))[0])[0] for f in interval_files]

        midpoints_list = []
        coverage_list = []

        for file in interval_files:
            midpoints, coverage = self.compute_binned_coverage(file, chromosome, bin_size, slide)
            midpoints_list.append(midpoints)
            coverage_list.append(coverage)
        self.plot_coverage(midpoints_list, coverage_list, labels, chromosome, output_dir)



class MutationDensityPlotter:
    def __init__(self, fai_file: str, verbose: bool = True):
        self.chrom_lengths = self._parse_fai(fai_file)
        self.verbose = verbose

    def _parse_fai(self, fai_file: str) -> dict:
        chrom_lengths = {}
        with open(fai_file) as f:
            for line in f:
                fields = line.strip().split('\t')
                chrom = fields[0]
                length = int(fields[1])
                chrom_lengths[chrom] = length
        return chrom_lengths

    def compute_mutation_density(
        self, mutation_file: str, chrom: str, bin_size: int, slide: int,
        mut_regex: Optional[re.Pattern] = None
    ) -> Tuple[List[int], List[int]]:
        chrom_length = self.chrom_lengths.get(chrom)
        df = pd.read_csv(mutation_file, compression='infer')
        df = df[df['chromosome'] == chrom]
        if mut_regex:
            df = df[df["mutation"].str.contains(mut_regex, regex=True, na=False)]

        starts = np.arange(0, chrom_length - bin_size + 1, slide)
        mutation_counts = [
            ((df['position'] >= start) & (df['position'] < start + bin_size)).sum()
            for start in starts
        ]
        midpoints = [start + bin_size // 2 for start in starts]
        return midpoints, mutation_counts

    def compute_coverage_for_normalization(
        self, interval_file: str, chrom: str, bin_size: int, slide: int
    ) -> List[float]:
        chrom_length = self.chrom_lengths.get(chrom)
        df = pd.read_csv(
            interval_file,
            sep='\t',
            compression='infer',
            dtype={"chromosome": str, "start": int, "end": int},
            header=0
        )
        df = df[df["chromosome"] == chrom].sort_values("start").reset_index(drop=True)

        starts = np.arange(0, chrom_length - bin_size + 1, slide)
        coverage = []

        current_index = 0
        n = len(df)

        for bin_start in starts:
            bin_end = bin_start + bin_size
            total_overlap = 0

            while current_index < n and df.at[current_index, "end"] <= bin_start:
                current_index += 1

            check_index = current_index
            while check_index < n and df.at[check_index, "start"] < bin_end:
                read_start = df.at[check_index, "start"]
                read_end = df.at[check_index, "end"]
                overlap_start = max(read_start, bin_start)
                overlap_end = min(read_end, bin_end)
                total_overlap += max(0, overlap_end - overlap_start)
                check_index += 1

            coverage.append(total_overlap / bin_size)

        return coverage

    def plot_mutation_density(
        self, midpoints_list: List[List[int]], values_list: List[List[float]],
        labels: List[str], chrom: str, output_path: str,
        normalized: bool = False, regex = ''
    ):
        plt.figure(figsize=(15, 5))
        for midpoints, values, label in zip(midpoints_list, values_list, labels):
            plt.plot(midpoints, values, label=label, lw=1.5)

        plt.title(f"{regex + ' ' if regex else ''}Mutation Density over {chrom}{' (normalized)' if normalized else ''}")
        plt.xlabel("Genomic Position")
        plt.ylabel("Mutation Density")
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(output_path, dpi=300)
        plt.close()

        log(f"Saved mutation density plot: {output_path}", self.verbose)


    def plot(
        self,
        mutation_dir: str,
        chromosome: str,
        output_dir: str,
        bin_size: int = 100000,
        slide: Optional[int] = None,
        coverage_dir: Optional[str] = None,
        mutation_category: Optional[str] = None
        ):
        if slide is None:
            slide = bin_size

        chrom = chromosome
        chrom_length = self.chrom_lengths.get(chrom)
        if chrom_length is None:
            raise ValueError(f"Chromosome {chrom} not found in FAI file.")

        mutation_files = sorted([
            os.path.join(mutation_dir, f)
            for f in os.listdir(mutation_dir)
            if f.endswith(".csv.gz")
        ])
        labels = [os.path.basename(f).replace("_mutations.csv.gz", "") for f in mutation_files]

        if coverage_dir:
            interval_files = sorted([
                os.path.join(coverage_dir, f)
                for f in os.listdir(coverage_dir)
                if f.endswith(".tsv") or f.endswith(".tsv.gz")
            ])
        else:
            interval_files = []

        mutation_regex = re.compile(mutation_category) if mutation_category else None

        midpoints_list = []
        values_list = []

        for i, mutation_file in enumerate(mutation_files):
            midpoints, mutations = self.compute_mutation_density(
                mutation_file, chrom, bin_size, slide, mutation_regex
            )

            if coverage_dir:
                coverage_file = interval_files[i]
                coverage = self.compute_coverage_for_normalization(
                    coverage_file, chrom, bin_size, slide
                )
                mutations = [m / c if c > 0 else 0 for m, c in zip(mutations, coverage)]

            midpoints_list.append(midpoints)
            values_list.append(mutations)

        plot_name = f"mutation_density_{chrom}"
        if mutation_category:
            plot_name += f"_{mutation_category}"
        if coverage_dir:
            plot_name += "_normalized"
        plot_name += ".png"

        output_path = os.path.join(output_dir, plot_name)
        self.plot_mutation_density(midpoints_list, values_list, labels, chrom, output_path, normalized=bool(coverage_dir), regex = mutation_category)

