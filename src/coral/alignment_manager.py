from collections import defaultdict
from io import TextIOWrapper
import os
import subprocess
import sys
import multiprocessing

import pysam
from .utils import run_cmd, log  
from typing import Optional, TextIO
import matplotlib.pyplot as plt


# Constants
LOW_MAPQ_THRESHOLD = 1
DEFAULT_MAPQ_THRESHOLD = 60
DEFAULT_OFFSET = 75


def log_to_file(log_path: Optional[str], message: str):
    if log_path:
        with open(log_path, 'a') as f:
            f.write(message + '\n')

def filter_sam(
    input_stream: TextIO,
    output_stream: TextIO,
    mapq_threshold: int = 60,
    mapq_hist_folder: Optional[str] = None,
    hist_name: str = None,
    verbose: bool = True,
    log_path: Optional[str] = None,
):
    total_reads = 0
    kept_reads = 0
    filtered_mapq = 0
    mapq_values = defaultdict(int)

    def write_summary():
        lines = [
            "Filter Summary:",
            f"  Total reads processed:   {total_reads}",
            f"  Reads kept:              {kept_reads}",
            f"  Filtered (low MAPQ):     {filtered_mapq}"
        ]
        if verbose:
            for line in lines:
                print(line)
        for line in lines:
            log_to_file(log_path, line)

    def plot_histogram():
        if mapq_hist_folder:
            scores = sorted(mapq_values.keys())
            counts = [mapq_values[score] for score in scores]

            plt.figure(figsize=(8, 5))
            plt.bar(scores, counts, color='steelblue', edgecolor='black', log=True)
            plt.title("MAPQ Score Distribution")
            plt.xlabel("MAPQ")
            plt.ylabel("Read Count (log scale)")
            plt.tight_layout()

            os.makedirs(mapq_hist_folder, exist_ok=True)
            out_path = os.path.join(mapq_hist_folder, hist_name)
            plt.savefig(out_path)
            msg = f"MAPQ histogram saved to {out_path}"
            if verbose:
                print(msg)
            log_to_file(log_path, msg)

    for i, line in enumerate(input_stream):
        if line.startswith('@'):
            output_stream.write(line)
            continue
        fields = line.split('\t')
        try:
            mapq = int(fields[4])
            mapq_values[mapq] += 1
            total_reads += 1
            if mapq >= mapq_threshold:
                output_stream.write(line)
                kept_reads += 1
            else:
                filtered_mapq += 1
        except Exception:
            msg = f"Invalid line {i+1}"
            print(msg, file=sys.stderr)
            log_to_file(log_path, msg)

    write_summary()
    plot_histogram()


def overlaps(read, other_reads):
    for other in other_reads:
        if read.reference_name != other.reference_name:
            continue
        if max(read.reference_start, other.reference_start) < min(read.reference_end, other.reference_end):
            return True
    return False

def with_continuity_filter_sam(
    input_stream,
    output_stream,
    low_mapq: int = 1,
    mapq_threshold: int = 60,
    mapq_hist_folder: Optional[str] = None,
    hist_name: str = None,
    verbose: bool = True,
    log_path: Optional[str] = None,
):
    total_reads = 0
    kept_reads = 0
    filtered_poor_mapping = 0
    filtered_mapq = 0
    filtered_disjoint = 0
    filtered_chrom = 0
    mapq_values = defaultdict(int)

    def write_summary():
        lines = [
            "Filter Summary:",
            f"  Total reads processed:   {total_reads}",
            f"  Reads kept:              {kept_reads}",
            f"  Filtered (poor mapping):     {filtered_poor_mapping}",
            f"  Filtered (low MAPQ):     {filtered_mapq}",
            f"  Filtered (disjoint):     {filtered_disjoint}",
            f"  Filtered (alt contigs):  {filtered_chrom}"
        ]
        if verbose:
            for line in lines:
                log(line, verbose)
        for line in lines:
            log_to_file(log_path, line)

    def plot_histogram():
        if mapq_hist_folder:
            scores = sorted(mapq_values.keys())
            counts = [mapq_values[score] for score in scores]
            plt.figure(figsize=(8, 5))
            plt.bar(scores, counts, color='steelblue', edgecolor='black', log=True)
            plt.title("MAPQ Score Distribution")
            plt.xlabel("MAPQ")
            plt.ylabel("Read Count (log scale)")
            plt.tight_layout()
            os.makedirs(mapq_hist_folder, exist_ok=True)
            out_path = os.path.join(mapq_hist_folder, hist_name)
            plt.savefig(out_path)
            msg = f"MAPQ histogram saved to {out_path}"
            log(msg, verbose)
            log_to_file(log_path, msg)

    prev_reads = []
    cur_reads = []
    next_reads = []
    next_read_name = None
    skip_contigs = {'Un', 'random', 'alt', 'fix', 'hap'}

    bamfile = pysam.AlignmentFile(input_stream, "rb")
    output_sam = pysam.AlignmentFile(output_stream, "wh", template=bamfile)

    for read in bamfile.fetch():
        total_reads += 1
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            filtered_poor_mapping += 1
            continue
        if read.mapping_quality < low_mapq:
            filtered_mapq += 1
            continue
        if next_read_name == read.query_name:
            next_reads.append(read)
        else:
            for cur_read in cur_reads:
                mapq_values[cur_read.mapping_quality] += 1
                chrom = cur_read.reference_name
                if any(keyword in chrom for keyword in skip_contigs):
                    filtered_chrom += 1
                    continue
                if cur_read.mapping_quality >= mapq_threshold or \
                   overlaps(cur_read, prev_reads) or overlaps(cur_read, next_reads):
                    output_sam.write(cur_read)
                    kept_reads += 1
                else:
                    filtered_disjoint += 1

            prev_reads, cur_reads, next_reads = cur_reads, next_reads, [read]
            next_read_name = read.query_name

    # Final group processing
    for cur_read in next_reads:
        total_reads += 1
        mapq_values[cur_read.mapping_quality] += 1
        chrom = cur_read.reference_name
        if any(keyword in chrom for keyword in skip_contigs):
            filtered_chrom += 1
            continue
        if cur_read.mapping_quality < low_mapq:
            filtered_mapq += 1
            continue
        if cur_read.mapping_quality >= mapq_threshold or overlaps(cur_read, cur_reads):
            output_sam.write(cur_read)
            kept_reads += 1
        else:
            filtered_disjoint += 1

    bamfile.close()
    output_sam.close()
    write_summary()
    plot_histogram()



class Aligner:
    def __init__(
        self,
        species_genome,
        reference_genome,
        base_output_dir,
        aligner_cmd=None,
        aligner_name=None,
        no_cache=False,
        cores = None,
        verbose=True
    ):
        self.species = species_genome.name
        self.reference = reference_genome.name
        self.output_dir = base_output_dir
        self.no_cache = no_cache
        self.verbose = verbose
        self.cores = cores if cores else multiprocessing.cpu_count()

        self.species_fasta = species_genome.fasta_path
        self.reference_fasta = reference_genome.fasta_path
        self.fastq = species_genome.fastq_path
        self.bam_dir = f"{self.output_dir}/BAMs"
        self.plots_dir = f"{self.output_dir}/Plots"
        self.raw_bam = f"{self.bam_dir}/{self.species}_to_{self.reference}_raw.bam"
        self.final_bam = f"{self.bam_dir}/{self.species}_to_{self.reference}.bam"
        self.hist_name = f"{self.species}_to_{self.reference}.png"
        self.log_path = self.final_bam.replace(".bam", ".log")

        os.makedirs(self.bam_dir, exist_ok=True)
        os.makedirs(self.plots_dir, exist_ok=True)

        if aligner_cmd:
            self.aligner_cmd_template = aligner_cmd
        elif aligner_name:
            self.aligner_cmd_template = self.get_aligner_cmd_from_name(aligner_name)
        else:
            log("No aligner specified. Defaulting to 'bwa'", self.verbose)
            self.aligner_cmd_template = self.get_aligner_cmd_from_name("bwa")

        self.validate_aligner_cmd()

    def get_aligner_cmd_from_name(self, name):
        commands = {
            "bwa": "bwa mem -t {cores} {ref} {fq}",
            "bwa-mem2": "bwa-mem2 mem -t {cores} {ref} {fq}",
            "minimap2": "minimap2 -t {cores} -ax sr {ref} {fq}",
            #"bbmap": "bbmap.sh ref={ref} threads={cores} in={fq} out=stdout.sam"
            "bbmap": "bbmap.sh ref={ref} threads={cores} in={fq} out=stdout.sam path={tmp}"

        }
        if name not in commands:
            raise ValueError(f"Unsupported aligner: {name}")
        return commands[name]

    def validate_aligner_cmd(self):
        required = ["{ref}", "{fq}", "{cores}"]
        for r in required:
            if r not in self.aligner_cmd_template:
                raise ValueError(f"--aligner-cmd must include placeholders: {', '.join(required)}")

    def align_streamed(self, mapq=60, low_mapq = 1, max_sort_mem=None, continuity = True):
        if os.path.exists(self.final_bam) and os.path.exists(self.final_bam + '.bai') and not self.no_cache:
            log(f"Streamed alignment already exists: {self.final_bam}", self.verbose)
            return self.final_bam

        cmd = self.aligner_cmd_template \
        .replace("{ref}", self.reference_fasta) \
        .replace("{fq}", self.fastq) \
        .replace("{cores}", str(self.cores)) \
        .replace("{tmp}", os.path.dirname(self.reference_fasta))

        log("Running full streaming alignment + filtering + sorting...", self.verbose)

        align_proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        sort_cmd = ["samtools", "sort", "-@", str(self.cores), "-o", self.final_bam]
        if max_sort_mem:
            sort_cmd += ["-m", str(max_sort_mem)]

        sort_proc = subprocess.Popen(sort_cmd, stdin=subprocess.PIPE)
        assert align_proc.stdout and sort_proc.stdin

        with TextIOWrapper(align_proc.stdout) as reader, TextIOWrapper(sort_proc.stdin, write_through=True, buffering=1) as writer:
            if continuity:
                with_continuity_filter_sam(
                input_stream=reader,
                output_stream=writer,
                low_mapq=low_mapq,
                mapq_threshold=mapq,
                mapq_hist_folder=self.plots_dir,
                hist_name=self.hist_name,
                verbose=self.verbose,
                log_path=self.log_path
                )
            else:
                filter_sam(
                    input_stream=reader,
                    output_stream=writer,
                    mapq_threshold=mapq,
                    mapq_hist_folder=self.plots_dir,
                    hist_name=self.hist_name,
                    verbose=self.verbose,
                    log_path=self.log_path
                )
        sort_proc.stdin.close()
        sort_proc.wait()
        run_cmd(["samtools", "index", self.final_bam])
        log(f"Finished (streamed): {self.final_bam}", self.verbose)
        return self.final_bam

    def align_disk_cached(self, mapq=60, low_mapq=1, continuity = True):
        # Step 1: Align and sort raw.bam
        if not os.path.exists(self.raw_bam) or self.no_cache:
            cmd = self.aligner_cmd_template \
            .replace("{ref}", self.reference_fasta) \
            .replace("{fq}", self.fastq) \
            .replace("{cores}", str(self.cores)) \
            .replace("{tmp}", os.path.dirname(self.reference_fasta))
            log("Running alignment to raw BAM...", self.verbose)

            tmp_raw_bam = self.raw_bam + ".tmp"
            run_cmd(f"{cmd} | samtools view -bS - > {tmp_raw_bam}", shell=True)
            os.rename(tmp_raw_bam, self.raw_bam) 

        else:
            log("Using cached raw BAM.", self.verbose)

        # Step 2: Filter from raw.bam to final.bam
        if not os.path.exists(self.final_bam) or self.no_cache:
            tmp_final_bam = self.final_bam + ".tmp"
            
            log("Filtering from raw BAM to final BAM...", self.verbose)
            view_proc = subprocess.Popen(["samtools", "view", "-h", self.raw_bam], stdout=subprocess.PIPE)
            sort_proc = subprocess.Popen(["samtools", "sort", "-@", str(self.cores), "-o", tmp_final_bam], stdin=subprocess.PIPE)
            assert view_proc.stdout and sort_proc.stdin

            if continuity:
                with_continuity_filter_sam(
                    input_stream=TextIOWrapper(view_proc.stdout),
                    output_stream=TextIOWrapper(sort_proc.stdin),
                    mapq_threshold=mapq,
                    low_mapq=low_mapq,
                    mapq_hist_folder=self.plots_dir,
                    hist_name=self.hist_name,
                    verbose=self.verbose,
                    log_path=self.log_path
                )
            else:
                filter_sam(
                    input_stream=TextIOWrapper(view_proc.stdout),
                    output_stream=TextIOWrapper(sort_proc.stdin),
                    mapq_threshold=mapq,
                    mapq_hist_folder=self.plots_dir,
                    hist_name=self.hist_name,
                    verbose=self.verbose,
                    log_path=self.log_path
                )
            
            sort_proc.stdin.close()
            sort_proc.wait()
            os.rename(tmp_final_bam, self.final_bam) 

        # Step 3: Index final.bam
        if (not os.path.exists(self.final_bam + '.bai') and not os.path.exists(self.final_bam + '.csi')) or self.no_cache:
            # run_cmd(["samtools", "index", self.final_bam])
            run_cmd(["samtools", "index", '-c', self.final_bam])

        log(f"Finished (disk-cached): {self.final_bam}", self.verbose)
        return self.final_bam


