import os
import subprocess
import shutil
from .utils import run_cmd, log
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

    def download(self):
        if os.path.exists(self.fasta_path) and not self.no_cache:
            log(f"Genome for {self.name} found in cache. Skipping download.", self.verbose)
            return

        os.makedirs(self.output_dir, exist_ok=True)
        temp_dir = os.path.join(os.path.dirname(self.output_dir), f"temp_{self.name}")
        os.makedirs(temp_dir, exist_ok=True)

        log(f"Downloading genome for {self.name} ({self.accession})", self.verbose)
        zip_path = os.path.join(temp_dir, f"{self.name}.zip")
        run_cmd([
            "datasets", "download", "genome", "accession", self.accession,
            "--filename", zip_path
        ], verbose=self.verbose)

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

        # Move the FASTA file to final location (efficient alternative to copy)
        shutil.move(fna_source, self.fasta_path)

        # Strip header and clean up
        log("Stripping FASTA headers (keeping only the first word)...", self.verbose)
        subprocess.run(["sed", "-i", "s/ .*$//", self.fasta_path], check=True)

        shutil.rmtree(temp_dir)



    def index(self, aligner="bwa"):
        log(f"Indexing genome for {self.name} with {aligner}...", self.verbose)

        aligner = (aligner or "").lower()

        # We only explicitly index for the BWA family.
        # minimap2 and bbmap build indices implicitly or in-memory.
        if aligner not in {"bwa", "bwa-mem2"}:
            # Still ensure FASTA index exists (used by samtools / pileup)
            fai = f"{self.fasta_path}.fai"
            if not os.path.exists(fai) or self.no_cache:
                if os.path.exists(fai):
                    os.remove(fai)
                run_cmd(["samtools", "faidx", self.fasta_path], verbose=self.verbose)
            log(f"No explicit index required for aligner '{aligner}'.", self.verbose)
            return

        # ---- BWA family explicit indexing ----

        if aligner == "bwa-mem2":
            # bwa-mem2 index files
            required_exts = ["bwt.2bit.64", "pac", "sa"]
            index_cmd = ["bwa-mem2", "index", self.fasta_path]
        else:
            # classic bwa index files
            required_exts = ["amb", "ann", "bwt", "pac", "sa"]
            index_cmd = ["bwa", "index", self.fasta_path]

        required_all = [f"{self.fasta_path}.{ext}" for ext in required_exts]
        required_all.append(f"{self.fasta_path}.fai")

        all_exist = all(os.path.exists(p) for p in required_all)

        if all_exist and not self.no_cache:
            log(f"Index files already exist for {self.name}. Skipping indexing.", self.verbose)
            return

        # Remove existing index files (safe for re-index)
        for ext in required_exts + ["fai"]:
            path = f"{self.fasta_path}.{ext}"
            if os.path.exists(path):
                os.remove(path)

        # Build aligner index
        run_cmd(index_cmd, verbose=self.verbose)

        # Always build samtools FASTA index
        run_cmd(["samtools", "faidx", self.fasta_path], verbose=self.verbose)

        log(f"Indexing complete for {self.name}", self.verbose)


    def generate_fragment_fastq(self, length=150, output_fastq = None, offset=75, force=False):
        output_fastq = output_fastq if output_fastq else os.path.join(self.output_dir, f"{self.name}.fastq")
        
        if os.path.exists(output_fastq) and not force:
            log(f"Fastq for {self.name} already exists. Skipping.", self.verbose)
            self.fastq_path = output_fastq
            return output_fastq


        def split_sequence(seq, offset, section_size):
            last_start = len(seq) - section_size
            sections = [(i, seq[i:i+section_size]) for i in range(0, last_start, offset)]
            if last_start > 0:
                sections.append((last_start, seq[last_start:]))
            return sections

        chromosomes = {
            record.id: record.seq
            for record in SeqIO.parse(self.fasta_path, "fasta")
        }

        with open(output_fastq, 'w') as out:
            for chrom_name, sequence in chromosomes.items():
                for i, (start, frag) in enumerate(split_sequence(sequence, offset, length)):
                    end = start + len(frag) - 1
                    out.write(f"@{chrom_name}_{start+1}_{end+1}\n")
                    out.write(f"{frag}\n+\n")
                    out.write(f"{'I' * len(frag)}\n")

        log(f"Wrote {output_fastq}", self.verbose)
        self.fastq_path = output_fastq
        return output_fastq

