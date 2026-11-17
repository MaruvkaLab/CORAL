from collections import defaultdict, deque
from io import TextIOWrapper
import os
import subprocess
import sys
import multiprocessing

import pysam
from .utils import run_cmd, log  
from typing import Optional, TextIO
import matplotlib.pyplot as plt


LOW_MAPQ_THRESHOLD = 1
DEFAULT_MAPQ_THRESHOLD = 60
DEFAULT_OFFSET = 75

# Built-in aligners: the command that indexes the reference (run once) and the alignment command.
# Placeholders: {ref} the reference FASTA, {refdir} its folder, {fq} the pseudo-reads, {cores} threads.
# minimap2's preset (-x sr) sets the index's k-mer and window sizes, so it is given in both commands.
ALIGNERS = {
    "bwa": {"index": "bwa index {ref}",
            "align": "bwa mem -t {cores} {ref} {fq}"},
    "bwa-mem2": {"index": "bwa-mem2 index {ref}",
                 "align": "bwa-mem2 mem -t {cores} {ref} {fq}"},
    "minimap2": {"index": "minimap2 -x sr -d {ref}.mmi {ref}",
                 "align": "minimap2 -t {cores} -ax sr {ref}.mmi {fq}"},
    "bbmap": {"index": "bbmap.sh ref={ref} path={refdir}",
              "align": "bbmap.sh path={refdir} threads={cores} in={fq} out=stdout.sam ordered=t"},
}


def aligner_commands(aligner_name="bwa", aligner_cmd=None, index_cmd=None):
    """(index command or None, alignment command): a built-in aligner, or a custom alignment command with its own
    optional index command. index_cmd also replaces a built-in aligner's index command."""
    if aligner_cmd:
        return index_cmd, aligner_cmd
    name = aligner_name or "bwa"
    if name not in ALIGNERS:
        raise ValueError(f"Unsupported aligner: {name} (built-in: {', '.join(ALIGNERS)}; or use --aligner-cmd)")
    return index_cmd or ALIGNERS[name]["index"], ALIGNERS[name]["align"]


def has_bam_index(bam_path):
    """True if the BAM has a .csi or .bai index.

    CORAL writes .csi (a .bai cannot address contigs over 512 Mb), but either one works."""
    return os.path.exists(bam_path + '.bai') or os.path.exists(bam_path + '.csi')


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
    """Keep reads with MAPQ >= mapq_threshold (no continuity rule)."""
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

    # Flush before the caller closes the pipe, or the last buffered reads are lost
    output_stream.flush()
    write_summary()
    plot_histogram()
    return {
        "total_reads": total_reads,
        "kept_reads": kept_reads,
        "filtered_low_mapq": filtered_mapq,
        "mapq_histogram": {str(q): n for q, n in sorted(mapq_values.items())},
    }


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
    continuity_run: int = 2,
):
    """Keep MAPQ >= mapq_threshold reads, and reads with low_mapq <= MAPQ < mapq_threshold
    that overlap on the reference with a neighbouring pseudo-read (upstream or downstream).

    Input must be in pseudo-read (query) order, as written by the aligner. continuity_run > 2
    requires a longer run of consecutive pseudo-reads, each overlapping the next.
    Secondary alignments and unplaced or alternative contigs are excluded."""
    if continuity_run < 2:
        raise ValueError("continuity_run must be at least 2")
    total_reads = 0
    kept_reads = 0
    kept_high_mapq = 0
    kept_rescued = 0
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
            f"  Kept (MAPQ >= threshold): {kept_high_mapq}",
            f"  Kept (continuity rescue): {kept_rescued}",
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

    # Each item is all alignments of one pseudo-read: before = already emitted, after = pending
    before = deque(maxlen=continuity_run - 1)
    after = deque()
    cur_reads = []
    cur_read_name = None
    skip_contigs = {'Un', 'random', 'alt', 'fix', 'hap'}

    bamfile = pysam.AlignmentFile(input_stream, "rb")
    try:
        output_sam = pysam.AlignmentFile(output_stream, "wh", template=bamfile)
    except Exception:
        bamfile.close()
        raise

    def run_length(read, fragments):
        # Longest chain of overlapping alignments from read through the following pseudo-reads
        chain = [[read]] + fragments
        run = [1] * len(chain[-1])
        for j in range(len(chain) - 2, -1, -1):
            run = [1 + max([n for other, n in zip(chain[j + 1], run) if overlaps(r, [other])], default=0)
                   for r in chain[j]]
        return run[0]

    def emit_first():
        nonlocal kept_reads, kept_high_mapq, kept_rescued, filtered_disjoint, filtered_chrom
        reads = after.popleft()
        for read in reads:
            mapq_values[read.mapping_quality] += 1
            chrom = read.reference_name
            if any(keyword in chrom for keyword in skip_contigs):
                filtered_chrom += 1
                continue
            if read.mapping_quality >= mapq_threshold:
                output_sam.write(read)
                kept_reads += 1
                kept_high_mapq += 1
                continue
            # Rescue if the run to the right, or both sides joined at this read, is long enough
            right = run_length(read, list(after))
            if right >= continuity_run or run_length(read, list(reversed(before))) + right - 1 >= continuity_run:
                output_sam.write(read)
                kept_reads += 1
                kept_rescued += 1
            else:
                filtered_disjoint += 1
        before.append(reads)

    try:
        for read in bamfile.fetch():
            total_reads += 1
            if read.is_unmapped or read.is_secondary:
                filtered_poor_mapping += 1
                continue
            if read.mapping_quality < low_mapq:
                filtered_mapq += 1
                continue
            if cur_read_name == read.query_name:
                cur_reads.append(read)
            else:
                if cur_reads:
                    after.append(cur_reads)
                    if len(after) == continuity_run:
                        emit_first()
                cur_reads = [read]
                cur_read_name = read.query_name

        if cur_reads:
            after.append(cur_reads)
        while after:
            emit_first()
    finally:
        bamfile.close()
        output_sam.close()

    write_summary()
    plot_histogram()
    return {
        "total_reads": total_reads,
        "kept_reads": kept_reads,
        "kept_high_mapq": kept_high_mapq,
        "kept_rescued": kept_rescued,
        "filtered_poor_mapping": filtered_poor_mapping,
        "filtered_low_mapq": filtered_mapq,
        "filtered_disjoint": filtered_disjoint,
        "filtered_alt_contig": filtered_chrom,
        "mapq_histogram": {str(q): n for q, n in sorted(mapq_values.items())},
    }



def annotate_sam(
    input_stream,
    output_stream,
    low_mapq: int = 1,
    tag: str = "ZN",
    verbose: bool = True,
    log_path: Optional[str] = None,
):
    """Write every primary mapped read with MAPQ >= low_mapq, tagged with the highest MAPQ
    of an overlapping neighbouring pseudo-read (-1 if none)."""
    total_reads = written = filtered_poor_mapping = filtered_mapq = filtered_chrom = 0
    skip_contigs = {'Un', 'random', 'alt', 'fix', 'hap'}

    bamfile = pysam.AlignmentFile(input_stream, "rb")
    try:
        output_sam = pysam.AlignmentFile(output_stream, "wh", template=bamfile)
    except Exception:
        bamfile.close()
        raise

    def best_overlap(read, others):
        best = -1
        for other in others:
            if other.reference_name == read.reference_name and \
               max(read.reference_start, other.reference_start) < min(read.reference_end, other.reference_end):
                best = max(best, other.mapping_quality)
        return best

    def emit(reads, before, after):
        nonlocal written, filtered_chrom
        for read in reads:
            if any(keyword in read.reference_name for keyword in skip_contigs):
                filtered_chrom += 1
                continue
            read.set_tag(tag, max(best_overlap(read, before), best_overlap(read, after)), value_type='i')
            output_sam.write(read)
            written += 1

    prev_reads, cur_reads, next_reads = [], [], []
    next_read_name = None
    try:
        for read in bamfile.fetch(until_eof=True):
            total_reads += 1
            if read.is_unmapped or read.is_secondary:
                filtered_poor_mapping += 1
                continue
            if read.mapping_quality < low_mapq:
                filtered_mapq += 1
                continue
            if next_read_name == read.query_name:
                next_reads.append(read)
            else:
                emit(cur_reads, prev_reads, next_reads)
                prev_reads, cur_reads, next_reads = cur_reads, next_reads, [read]
                next_read_name = read.query_name
        emit(cur_reads, prev_reads, next_reads)
        emit(next_reads, cur_reads, [])
    finally:
        bamfile.close()
        output_sam.close()

    stats = {
        "total_reads": total_reads,
        "written_reads": written,
        "filtered_poor_mapping": filtered_poor_mapping,
        "filtered_low_mapq": filtered_mapq,
        "filtered_alt_contig": filtered_chrom,
    }
    for line in ["Annotate summary:"] + [f"  {k}: {v}" for k, v in stats.items()]:
        log(line, verbose)
        log_to_file(log_path, line)
    return stats


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
        verbose=True,
        annotate=False,
        index_cmd=None
    ):
        self.species = species_genome.name
        self.reference = reference_genome.name
        self.output_dir = base_output_dir
        self.no_cache = no_cache
        self.verbose = verbose
        self.cores = cores if cores else multiprocessing.cpu_count()

        self.species_fasta = species_genome.fasta_path
        self.reference_genome = reference_genome
        self.reference_fasta = reference_genome.fasta_path
        self.fastq = species_genome.fastq_path
        self.bam_dir = f"{self.output_dir}/BAMs"
        self.plots_dir = f"{self.output_dir}/Plots"
        self.raw_bam = f"{self.bam_dir}/{self.species}_to_{self.reference}_raw.bam"
        # annotate: annotate_sam replaces the MAPQ filter; own BAM name so the caches never mix
        self.annotate = annotate
        suffix = "_annotated" if annotate else ""
        self.final_bam = f"{self.bam_dir}/{self.species}_to_{self.reference}{suffix}.bam"
        self.hist_name = f"{self.species}_to_{self.reference}.png"
        self.log_path = self.final_bam.replace(".bam", ".log")
        self.filter_stats = None  # populated when this run actually filters (not when cached)

        os.makedirs(self.bam_dir, exist_ok=True)
        os.makedirs(self.plots_dir, exist_ok=True)

        self.index_cmd, self.aligner_cmd_template = aligner_commands(aligner_name, aligner_cmd, index_cmd)
        self.validate_aligner_cmd()

    def validate_aligner_cmd(self):
        cmd = self.aligner_cmd_template
        if "{fq}" not in cmd or "{cores}" not in cmd or not ("{ref}" in cmd or "{refdir}" in cmd):
            raise ValueError("--aligner-cmd must include {fq}, {cores}, and {ref} or {refdir}")

    def _align_command(self):
        """The alignment command, after making sure the reference is indexed (a no-op once it is)."""
        self.reference_genome.index(self.index_cmd, self.cores)
        ref = self.reference_fasta
        return (self.aligner_cmd_template.replace("{refdir}", os.path.dirname(os.path.abspath(ref)))
                .replace("{ref}", ref).replace("{fq}", self.fastq).replace("{cores}", str(self.cores)))

    def align_streamed(self, mapq=60, low_mapq = 1, max_sort_mem=None, continuity = True, continuity_run=2):
        if os.path.exists(self.final_bam) and has_bam_index(self.final_bam) and not self.no_cache:
            log(f"Streamed alignment already exists: {self.final_bam}", self.verbose)
            return self.final_bam

        cmd = self._align_command()

        log("Running full streaming alignment + filtering + sorting...", self.verbose)

        align_proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        sort_cmd = ["samtools", "sort", "-@", str(self.cores), "-o", self.final_bam]
        if max_sort_mem:
            sort_cmd += ["-m", str(max_sort_mem)]

        sort_proc = subprocess.Popen(sort_cmd, stdin=subprocess.PIPE)
        assert align_proc.stdout and sort_proc.stdin

        with TextIOWrapper(align_proc.stdout) as reader, TextIOWrapper(sort_proc.stdin, write_through=True) as writer:
            if self.annotate:
                self.filter_stats = annotate_sam(
                    input_stream=reader,
                    output_stream=writer,
                    low_mapq=low_mapq,
                    verbose=self.verbose,
                    log_path=self.log_path
                )
            elif continuity:
                self.filter_stats = with_continuity_filter_sam(
                input_stream=reader,
                output_stream=writer,
                low_mapq=low_mapq,
                mapq_threshold=mapq,
                mapq_hist_folder=self.plots_dir,
                hist_name=self.hist_name,
                verbose=self.verbose,
                log_path=self.log_path,
                continuity_run=continuity_run
                )
            else:
                self.filter_stats = filter_sam(
                    input_stream=reader,
                    output_stream=writer,
                    mapq_threshold=mapq,
                    mapq_hist_folder=self.plots_dir,
                    hist_name=self.hist_name,
                    verbose=self.verbose,
                    log_path=self.log_path
                )
        sort_proc.stdin.close()
        if align_proc.wait() != 0 or sort_proc.wait() != 0:
            raise RuntimeError(f"Streamed alignment failed (aligner {align_proc.returncode}, samtools sort {sort_proc.returncode})")
        run_cmd(["samtools", "index", '-c', self.final_bam])
        log(f"Finished (streamed): {self.final_bam}", self.verbose)
        return self.final_bam

    def align_disk_cached(self, mapq=60, low_mapq=1, continuity = True, continuity_run=2):
        # Step 1: align to an unsorted raw BAM (query order, as the filters need)
        if not os.path.exists(self.raw_bam) or self.no_cache:
            cmd = self._align_command()
            log("Running alignment to raw BAM...", self.verbose)

            tmp_raw_bam = self.raw_bam + ".tmp"
            # pipefail: an aligner that dies after writing the header must fail the run, not leave a truncated BAM
            run_cmd(["bash", "-o", "pipefail", "-c", f"{cmd} | samtools view -bS - > {tmp_raw_bam}"])
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

            if self.annotate:
                self.filter_stats = annotate_sam(
                    input_stream=TextIOWrapper(view_proc.stdout),
                    output_stream=TextIOWrapper(sort_proc.stdin),
                    low_mapq=low_mapq,
                    verbose=self.verbose,
                    log_path=self.log_path
                )
            elif continuity:
                self.filter_stats = with_continuity_filter_sam(
                    input_stream=TextIOWrapper(view_proc.stdout),
                    output_stream=TextIOWrapper(sort_proc.stdin),
                    mapq_threshold=mapq,
                    low_mapq=low_mapq,
                    mapq_hist_folder=self.plots_dir,
                    hist_name=self.hist_name,
                    verbose=self.verbose,
                    log_path=self.log_path,
                    continuity_run=continuity_run
                )
            else:
                self.filter_stats = filter_sam(
                    input_stream=TextIOWrapper(view_proc.stdout),
                    output_stream=TextIOWrapper(sort_proc.stdin),
                    mapq_threshold=mapq,
                    mapq_hist_folder=self.plots_dir,
                    hist_name=self.hist_name,
                    verbose=self.verbose,
                    log_path=self.log_path
                )
            
            sort_proc.stdin.close()
            if view_proc.wait() != 0 or sort_proc.wait() != 0:
                raise RuntimeError(f"Filtering {self.raw_bam} failed (samtools view {view_proc.returncode}, samtools sort {sort_proc.returncode})")
            os.rename(tmp_final_bam, self.final_bam) 

        # Step 3: Index final.bam
        if not has_bam_index(self.final_bam) or self.no_cache:
            run_cmd(["samtools", "index", '-c', self.final_bam])

        log(f"Finished (disk-cached): {self.final_bam}", self.verbose)
        return self.final_bam


