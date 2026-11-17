import os
import sys
import subprocess
import shutil
from itertools import chain
from .utils import run_cmd, run_cmd_retry, log
from Bio import SeqIO




class Genome:
    def __init__(self, name, accession, output_dir, fasta_path = None, no_cache=False, verbose = True):
        self.name = name
        self.accession = accession
        self.output_dir = os.path.join(output_dir, name)
        self.no_cache = no_cache
        self.fasta_path = fasta_path if fasta_path else os.path.join(self.output_dir, f"{name}.fasta")
        self.fastq_path = None
        self.verbose = verbose
        self.n_contigs = None   # set when the FASTA is parsed for fragmentation
        self.total_bp = None
        self.masked_bp = self.fragments_written = self.fragments_skipped_repeat = None

    def download(self):
        if os.path.exists(self.fasta_path) and not self.no_cache:
            log(f"Genome for {self.name} found in cache. Skipping download.", self.verbose)
            return

        os.makedirs(self.output_dir, exist_ok=True)
        temp_dir = os.path.join(os.path.dirname(self.output_dir), f"temp_{self.name}")
        os.makedirs(temp_dir, exist_ok=True)

        log(f"Downloading genome for {self.name} ({self.accession})", self.verbose)
        zip_path = os.path.join(temp_dir, f"{self.name}.zip")
        download_cmd = [
            "datasets", "download", "genome", "accession", self.accession,
            "--filename", zip_path
        ]
        # No progress bar when stderr is not a TTY (e.g., cluster jobs)
        if not sys.stderr.isatty():
            download_cmd.append("--no-progressbar")
        # Raises NCBI's per-IP rate limit (3 -> 10 req/s) when a key is available
        api_key = os.environ.get("NCBI_API_KEY")
        if api_key:
            download_cmd += ["--api-key", api_key]
        run_cmd_retry(download_cmd, verbose=self.verbose)

        log("Unzipping genome...", self.verbose)
        run_cmd(["unzip", "-q", "-o", zip_path, "-d", temp_dir], verbose=self.verbose)


        fna_source = None
        for root, _, files in os.walk(temp_dir):
            for file in files:
                if file.endswith(".fna"):
                    fna_source = os.path.join(root, file)
                    break

        if not fna_source:
            log("Error: No FASTA (.fna) file found in the downloaded data.", self.verbose)
            shutil.rmtree(temp_dir)
            raise FileNotFoundError("Missing FASTA")

        shutil.move(fna_source, self.fasta_path)

        log("Stripping FASTA headers (keeping only the first word)...", self.verbose)
        subprocess.run(["sed", "-i", "s/ .*$//", self.fasta_path], check=True)

        shutil.rmtree(temp_dir)



    def index(self, index_cmd=None, cores=1, force=False):
        """Index this genome as a reference: the samtools .fai and, with index_cmd (see alignment_manager.ALIGNERS),
        the aligner's index, each built once. <fasta>.coral_index records the command that built the index, so an
        index from another command, or an unfinished one, is rebuilt."""
        if force or not os.path.exists(f"{self.fasta_path}.fai"):
            run_cmd(["samtools", "faidx", self.fasta_path], verbose=self.verbose)
        if not index_cmd:
            return
        marker = f"{self.fasta_path}.coral_index"
        if not force and os.path.exists(marker) and open(marker).read() == index_cmd:
            log(f"Index already built for {self.name}.", self.verbose)
            return
        if os.path.exists(marker):
            os.remove(marker)
        cmd = (index_cmd.replace("{refdir}", os.path.dirname(os.path.abspath(self.fasta_path)))
               .replace("{ref}", self.fasta_path).replace("{cores}", str(cores)))
        run_cmd(["bash", "-o", "pipefail", "-c", cmd], verbose=self.verbose)
        with open(marker, "w") as f:
            f.write(index_cmd)

    def generate_fragment_fastq(self, length=150, output_fastq = None, offset=75, force=False,
                                repeat_mask=None, mask_frac=0.5):
        """Fragment the genome into pseudo-reads (150 bp at a 75 bp step by default) and write a FASTQ."""
        output_fastq = output_fastq if output_fastq else os.path.join(self.output_dir, f"{self.name}.fastq")

        if os.path.exists(output_fastq) and not force:
            log(f"Fastq for {self.name} already exists. Skipping.", self.verbose)
            self.fastq_path = output_fastq
            return output_fastq


        # With a repeat mask, skip fragments that are >= mask_frac masked (they mostly
        # mis-map to paralogs). Read names carry the 1-based inclusive span.
        self.n_contigs = self.total_bp = 0
        written = skipped = 0
        with open(output_fastq + ".tmp", 'w') as out:
            for record in SeqIO.parse(self.fasta_path, "fasta"):
                chrom_name = record.id
                sequence = str(record.seq)
                self.n_contigs += 1
                self.total_bp += len(sequence)

                last_start = len(sequence) - length
                starts = range(0, last_start, offset)
                if last_start >= 0:   # last fragment ends at the contig end; shorter contigs give none
                    starts = chain(starts, (last_start,))

                for start in starts:
                    frag = sequence[start:start + length]
                    end = start + len(frag) - 1
                    if repeat_mask and repeat_mask.masked_fraction(chrom_name, start + 1, end + 1) >= mask_frac:
                        skipped += 1
                        continue
                    out.write(f"@{chrom_name}_{start+1}_{end+1}\n")
                    out.write(f"{frag}\n+\n")
                    out.write(f"{'I' * len(frag)}\n")
                    written += 1

        os.replace(output_fastq + ".tmp", output_fastq)
        self.fragments_written, self.fragments_skipped_repeat = written, skipped
        self.masked_bp = repeat_mask.n_masked_bases if repeat_mask else None
        if repeat_mask:
            log(f"Wrote {output_fastq} ({written} fragments; {skipped} skipped as "
                f">={mask_frac:.0%} repeat)", self.verbose)
        else:
            log(f"Wrote {output_fastq}", self.verbose)
        self.fastq_path = output_fastq
        return output_fastq

