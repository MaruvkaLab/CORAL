# species_mutation_extraction/mutextractor/pipeline.py

import json
import os
import re
import subprocess
import time
import gc

try:
    import resource
except ImportError:
    resource = None

import pandas as pd
from .cleanup_manager import PipelineCleaner
from .genome_manager import Genome
from .alignment_manager import Aligner
from .multiple_species_mutation_extractor_manager import MultipleSpeciesMutationExtractor
from .mutation_extractor_manager import FiveMerExtractor, MutationExtractor, MutationNormalizer, ParallelMutationExtractor, TripletExtractor
from .pileup_manager import Pileup
from .plot_utils import CoveragePlotter, MutationDensityPlotter, MutationSpectraPlotter
from .utils import get_top_n_chromosomes, log
from .repeat_masker import build_mask
import psutil
import pysam


def _cpu_seconds():
    """CPU seconds used by this process and all its children, or None."""
    if resource is None:
        return None
    s = resource.getrusage(resource.RUSAGE_SELF)
    c = resource.getrusage(resource.RUSAGE_CHILDREN)
    return s.ru_utime + s.ru_stime + c.ru_utime + c.ru_stime


def _tool_versions():
    def pick(name, lines):
        if not lines:
            return None
        if name == "bwa-mem2":  # bare version on stdout, loader chatter on stderr
            return next((l for l in lines if re.match(r'^\d+\.\d+', l)), lines[0])
        # bwa prints "Version: x"; samtools prints "samtools x" first
        return next((l for l in lines if "ersion" in l), lines[0])

    versions = {}
    for name, cmd in (("samtools", ["samtools", "--version"]),
                      ("bwa", ["bwa"]),
                      ("bwa-mem2", ["bwa-mem2", "version"])):
        try:
            r = subprocess.run(cmd, capture_output=True, text=True)
            lines = [l.strip() for l in ((r.stdout or "") + (r.stderr or "")).splitlines() if l.strip()]
            found = pick(name, lines)
            if found:
                versions[name] = found
        except Exception:
            pass
    try:
        from importlib.metadata import version
        versions["coral"] = version("coral")
    except Exception:
        pass
    return versions


class MutationExtractionPipeline:
    def __init__(self, 
                 species_list,
                 outgroup,
                 aligner_name="bwa-mem2", 
                 aligner_cmd=None,
                 base_output_dir="../Output", 
                 no_cache = False,
                 verbose = True, 
                 run_id = None,
                 **kwargs):
        if len(species_list) != 2:
            raise ValueError(
                f"MutationExtractionPipeline needs exactly 2 ingroup species, got {len(species_list)}. "
                "Use MultiSpeciesMutationPipeline for more.")
        self.species_list = species_list  # list of (name, accession)
        self.outgroup = outgroup          # (name, accession)
        self.aligner_name = aligner_name
        self.aligner_cmd = aligner_cmd
        self.run_id = run_id
        if run_id is None:
            self.run_id = '__'.join([species[0] for species in [outgroup] + species_list])
        s = kwargs.get("suffix")
        if s:
            suffix = "_" + str(s)
        else:
            suffix = ""
        self.output_dir = f"{base_output_dir}/{self.run_id}{suffix}"
        self.params = kwargs

        # Will hold references to internal data
        self.reference = None
        self.genomes = []
        self.alignments = []
        self.align_times = {}
        self.genome_stats = {}
        self.reference_mask = None    # RepeatMask on the outgroup, if --repeat-mask
        self.verbose = verbose
        self.no_cache = no_cache

    def _build_repeat_mask(self, genome):
        """Build (cached) a repeat mask for `genome` if --repeat-mask is on, else
        None. Uses WindowMasker by default (library-free); RepeatMasker if selected."""
        if not self.params.get("repeat_mask", False):
            return None
        return build_mask(
            genome.fasta_path, genome.output_dir,
            tool=self.params.get("repeat_masker", "windowmasker"),
            species=self.params.get("repeat_species"),
            cores=self.params.get("cores", 1) or 1,
            dust=self.params.get("repeat_dust", True),
            no_cache=self.no_cache, verbose=self.verbose)

    
    def run(self):
        log("Starting mutation extraction pipeline...", self.verbose)
        timings = {}
        cpu_times = {}
        memory_log = {}
        process = psutil.Process(os.getpid())
        start_pipeline = time.time()

        def get_memory():
            return round(process.memory_info().rss / (1024 ** 2), 2)  # In MB

        def timed_stage(stage_name, func):
            log(f"--- Starting: {stage_name} ---", self.verbose)
            mem_before = get_memory()
            start = time.time()
            cpu_before = _cpu_seconds()
            func()
            gc.collect()  # Clean up memory after each stage
            end = time.time()
            mem_after = get_memory()
            cpu_after = _cpu_seconds()

            timings[stage_name] = round(end - start, 2)
            if cpu_before is not None:
                cpu_times[stage_name] = round(cpu_after - cpu_before, 2)
            memory_log[stage_name] = {"start_MB": mem_before, "end_MB": mem_after}
            log(f"{stage_name} completed in {timings[stage_name]} seconds", self.verbose)
            log(f"Memory usage: {mem_before} → {mem_after} MB", self.verbose)

        timed_stage("Download and Fragment Genomes", self.download_index_and_fragment_genomes)
        timed_stage("Align Species", self.align_species)
        timed_stage("Generate Pileup", self.generate_pileup)
        timed_stage("Extract Mutations and Triplets", self.extract_mutations_and_triplets)
        timed_stage("Extract Intervals", self.extract_intervals)
        if self.params.get("plots", True):
            timed_stage("Run Plots", self.run_plots)
        self.genome_stats = self._collect_genome_stats()  # before cleanup removes the FASTAs
        timed_stage("Cleanup files", self.cleanup)

        total_runtime = round(time.time() - start_pipeline, 2)
        timings["Total Runtime"] = total_runtime
        memory_log["Total Runtime"] = {"final_MB": get_memory()}

        timing_path = os.path.join(self.output_dir, "pipeline_timings.json")
        with open(timing_path, "w") as f:
            json.dump({"timings": timings, "memory": memory_log}, f, indent=2)

        try:
            self._write_run_summary(timings, cpu_times)
        except Exception as e:  # diagnostics must never fail a completed run
            log(f"Warning: could not write run summary: {e}", self.verbose)

        log(f"Timing and memory info saved to: {timing_path}", self.verbose)
        log("Pipeline completed successfully.", self.verbose)


    def _collect_genome_stats(self):
        stats = {}
        for genome in ([self.reference] if self.reference else []) + self.genomes:
            entry = {"accession": genome.accession}
            fai = genome.fasta_path + ".fai"
            if os.path.exists(fai):
                with open(fai) as f:
                    lengths = [int(line.split('\t')[1]) for line in f if line.strip()]
                entry["contigs"] = len(lengths)
                entry["total_bp"] = sum(lengths)
            elif genome.total_bp is not None:
                entry["contigs"] = genome.n_contigs
                entry["total_bp"] = genome.total_bp
            elif os.path.exists(genome.fasta_path):
                entry["fasta_bytes"] = os.path.getsize(genome.fasta_path)
            stats[genome.name] = entry
        return stats

    def _write_run_summary(self, timings, cpu_times):
        path = os.path.join(self.output_dir, "run_summary.json")
        data = {}
        if os.path.exists(path):
            try:
                with open(path) as f:
                    data = json.load(f)
            except (ValueError, OSError):
                data = {}

        data["run"] = {
            "run_id": self.run_id,
            "outgroup": self.outgroup[0],
            "taxa": [name for name, _ in self.species_list],
            "parameters": {k: v for k, v in self.params.items()
                           if v is None or isinstance(v, (str, int, float, bool))},
            "versions": _tool_versions(),
        }
        data["genomes"] = getattr(self, "genome_stats", {})
        alignment = {a.species: a.filter_stats for a in self.alignments if a.filter_stats}
        if alignment:
            data["alignment"] = alignment
        data["timings"] = {
            "wall_seconds": timings,
            "cpu_seconds": cpu_times,
            "alignment_wall_seconds": getattr(self, "align_times", {}),
        }

        tmp = path + ".tmp"
        with open(tmp, "w") as f:
            json.dump(data, f, indent=2)
        os.replace(tmp, path)
        log(f"Run summary saved to: {path}", self.verbose)

    def download_index_and_fragment_genomes(self):
        # log("Downloading, indexing, and fragmenting genomes...", self.verbose)

        # Reference genome (outgroup)
        ref_name, ref_acc = self.outgroup
        self.reference = Genome(
            name=ref_name,
            accession=ref_acc,
            output_dir=self.output_dir,
            no_cache=self.no_cache,
            verbose=self.verbose
        )
        self.reference.download()
        self.reference.index(aligner=self.aligner_name)
        # Outgroup repeat mask -> used to drop calls at reference repeat positions.
        self.reference_mask = self._build_repeat_mask(self.reference)

        # Ingroup genomes
        for name, acc in self.species_list:
            genome = Genome(
                name=name,
                accession=acc,
                output_dir=self.output_dir,
                no_cache=self.no_cache,
                verbose=self.verbose
            )
            genome.download()
            # Species repeat mask -> don't even generate pseudo-reads from repeats.
            species_mask = self._build_repeat_mask(genome)
            genome.generate_fragment_fastq(
                length=self.params.get("fragment_length", 150),
                offset=self.params.get("fragment_offset", 75),
                force=self.no_cache,
                repeat_mask=species_mask,
                mask_frac=self.params.get("repeat_mask_frac", 0.5))
            self.genomes.append(genome)

    def align_species(self):
        # log("Aligning species to reference...", self.verbose)

        for genome in self.genomes:
            aligner = Aligner(
                species_genome=genome,
                reference_genome=self.reference,
                base_output_dir=self.output_dir,
                aligner_cmd=self.aligner_cmd,
                aligner_name=self.aligner_name,              
                cores=self.params.get("cores", None),
                verbose=self.verbose,
            )

            t0 = time.time()
            if self.params.get("streamed", False):
                aligner.align_streamed(
                    mapq=self.params.get("mapq", 60),
                    low_mapq=self.params.get("low_mapq", 1),
                    max_sort_mem=self.params.get("max_samtools_mem", None),
                    continuity=self.params.get("continuity", True)
                )
            else:
                aligner.align_disk_cached(
                    mapq=self.params.get("mapq", 60),
                    low_mapq=self.params.get("low_mapq", 1),
                    continuity=self.params.get("continuity", True)
                )
            self.align_times[genome.name] = round(time.time() - t0, 2)
            self.alignments.append(aligner)
            


    def generate_pileup(self):
        # log("Generating pileup from alignments...", self.verbose)

        # Ensure aligners were run and final BAMs are available
        for aligner in self.alignments:
            if not os.path.exists(aligner.final_bam):
                raise FileNotFoundError(f"BAM not found: {aligner.final_bam}")

        pileup_generator = Pileup(
            outgroup=self.reference,
            aligners=self.alignments,
            base_output_dir=self.output_dir,
            run_id=self.run_id,
            no_cache=self.no_cache,
            verbose=self.verbose
        )
        
        self.pileup = pileup_generator
        self.pileup_path = pileup_generator.pileup_path

        cores = self.params.get("cores")
        parallel_extract = bool(cores) and cores > 1
        # In parallel mode the extractor generates per-chromosome pileups directly
        # from the indexed BAMs, so the whole-genome pileup is only needed for the
        # serial scan or the (opt-in) 5-mer pass.
        if not parallel_extract or self.params.get("five_mer", False):
            self.pileup_path = pileup_generator.generate()
        else:
            log("Parallel extraction enabled: skipping whole-genome pileup "
                "(per-chromosome pileups are generated during extraction).", self.verbose)


    def extract_mutations_and_triplets(self):
        # log("Extracting 3mer mutations and triplets from pileup...", self.verbose)
        cores = self.params.get("cores")
        mut_dir = os.path.join(self.output_dir, 'Mutations')
        trip_dir = os.path.join(self.output_dir, 'Triplets')
        # cores>1 -> chromosome-parallel extraction via per-chromosome mpileup
        # (byte-identical output); otherwise the unchanged serial scan over the
        # whole-genome pileup (no added overhead on a single core).
        if cores and cores > 1:
            mutation_extractor = ParallelMutationExtractor(
                reference=self.reference.name,
                taxon1=self.genomes[0].name,
                taxon2=self.genomes[1].name,
                ref_fasta=self.reference.fasta_path,
                bams=[aligner.final_bam for aligner in self.alignments],
                mutation_output_dir=mut_dir,
                triplet_output_dir=trip_dir,
                fai_path=self.reference.fasta_path + ".fai",
                cores=cores,
                no_full_mutations=False,
                no_cache=self.no_cache,
                verbose=self.verbose,
                ref_mask=self.reference_mask)
        else:
            mutation_extractor = MutationExtractor(
                reference=self.reference.name,
                taxon1=self.genomes[0].name,
                taxon2=self.genomes[1].name,
                pileup_file=self.pileup_path,
                mutation_output_dir=mut_dir,
                triplet_output_dir=trip_dir,
                no_full_mutations=False,
                no_cache=self.no_cache,
                verbose=self.verbose,
                ref_mask=self.reference_mask)
        mutation_extractor.extract()

        # 5-mer extraction is an extra full pass over the pileup and is not used
        # by the standard (96-category, trinucleotide) outputs or normalization,
        # so it is opt-in. Enable with five_mer=True (CLI: --five-mer).
        if self.params.get("five_mer", False):
            fivemer_extractor = FiveMerExtractor(reference=self.reference.name,
                                  taxon1=self.genomes[0].name,
                                  taxon2=self.genomes[1].name,
                                  pileup_file=self.pileup_path,
                                  output_dir=os.path.join(self.output_dir, 'Mutations'),
                                  no_cache=self.no_cache,
                                  verbose=self.verbose)
            fivemer_extractor.extract()

        # log("Extracting triplets from pileup...", self.verbose)
        # triplet_extractor = TripletExtractor(reference=self.reference.name,
        #                       taxon1=self.genomes[0].name,
        #                       taxon2=self.genomes[1].name,
        #                       pileup_file=self.pileup_path,
        #                       output_dir=os.path.join(self.output_dir, 'Triplets'),
        #                       no_cache=False,
        #                       verbose=self.verbose)
        # triplet_extractor.extract()

        normalizer = MutationNormalizer(
            input_dir=self.output_dir,
            output_dir= os.path.join(self.output_dir, "Tables"),
            verbose=self.verbose,
            divergence_time= self.params.get("divergence_time", None),
        )
        normalizer.normalize()

    def _extract_bam_intervals(self, input_bam, output_dir, assume_sorted=False, merge=False, no_cache=False):
            os.makedirs(output_dir, exist_ok=True)

            base_name = os.path.basename(input_bam).rsplit(".", 1)[0]
            output_file = os.path.join(output_dir, f"{base_name}_intervals.tsv.gz")

            if os.path.exists(output_file) and not no_cache:
                log(f"Intervals already exist: {output_file}", self.verbose)
                return output_file

            def extract_raw_intervals(bamfile):
                intervals = []
                for read in bamfile.fetch():
                    if not read.is_unmapped:
                        chrom = bamfile.get_reference_name(read.reference_id)
                        intervals.append((chrom, read.reference_start, read.reference_end))
                return intervals

            def extract_intervals_sorted(bamfile):
                merged = []
                for read in bamfile.fetch():
                    if read.is_unmapped:
                        continue
                    chrom = bamfile.get_reference_name(read.reference_id)
                    start = read.reference_start
                    end = read.reference_end
                    if merged and merged[-1][0] == chrom and merged[-1][2] >= start:
                        merged[-1] = (chrom, merged[-1][1], max(merged[-1][2], end))
                    else:
                        merged.append((chrom, start, end))
                return merged

            def merge_intervals(intervals):
                merged = []
                for chrom, start, end in sorted(intervals):
                    if merged and merged[-1][0] == chrom and merged[-1][2] >= start:
                        merged[-1] = (chrom, merged[-1][1], max(merged[-1][2], end))
                    else:
                        merged.append((chrom, start, end))
                return merged

            with pysam.AlignmentFile(input_bam, "rb") as bamfile:
                intervals = (
                    extract_intervals_sorted(bamfile)
                    if assume_sorted
                    else merge_intervals(extract_raw_intervals(bamfile)) if merge
                    else extract_raw_intervals(bamfile)
                )

            df = pd.DataFrame(intervals, columns=["chromosome", "start", "end"])
            df.to_csv(output_file, sep='\t', index=False, compression="gzip")

            log(f"Intervals written to: {output_file}", self.verbose)
            return output_file
    
    def extract_intervals(self):
        for bam in self.alignments:
            self._extract_bam_intervals(bam.final_bam, os.path.join(self.output_dir, 'Intervals'),
                                        no_cache=self.no_cache)

    def run_plots(self):
        spectra_plotter = MutationSpectraPlotter()
        spectra_plotter.plot(tables_dir = os.path.join(self.output_dir, 'Tables'))
        fai_file = self.reference.fasta_path + '.fai'
        coverage_plotter = CoveragePlotter(fai_file=fai_file)
        mutation_density_plotter = MutationDensityPlotter(fai_file=fai_file)

        top_chroms = get_top_n_chromosomes(fai_file, n=3)
        log("Plotting coverage and mutation density for top chromosomes...", self.verbose)
        for chrom in top_chroms:
            log(f"Plotting for {chrom}...", self.verbose)

            coverage_plotter.plot(interval_dir=os.path.join(self.output_dir, 'Intervals'),
                                 chromosome=chrom,
                                 output_dir=os.path.join(self.output_dir, 'Plots', f"coverage_{chrom}.png"))

            mutation_density_plotter.plot(mutation_dir=os.path.join(self.output_dir, 'Mutations'),
                                 chromosome=chrom,
                                 output_dir=os.path.join(self.output_dir, 'Plots'))
            
            mutation_density_plotter.plot(mutation_dir=os.path.join(self.output_dir, 'Mutations'),
                                 chromosome=chrom,
                                 output_dir=os.path.join(self.output_dir, 'Plots'),
                                 mutation_category = r"[ACTG][C>T]G")
            
    def cleanup(self):
        cleaner = PipelineCleaner(self.genomes + [self.reference], self.alignments, self.pileup, base_dir=self.output_dir, verbose=self.verbose)
        cleaner.run(bams=True, pileup=True, genomes=True)


from .multiple_species_utils import (
    annotate_list_with_indices,
    parse_species_accession_from_newick,
    annotate_tree_with_indices,
    save_annotated_tree,
)
from . import parallel
from .run_phylip import run_phylip, check_phylip_available


class MultiSpeciesMutationPipeline:
    def __init__(
        self,
        newick_tree = None,
        species_list = None,
        base_output_dir="../Output",
        run_id=None,
        outgroup=None,
        aligner_name="bwa-mem2", 
        aligner_cmd=None,
        no_cache=False,
        verbose=True,
        cores=None,
        no_parallel=False,
        align_jobs=None,
        scan_jobs=None,
        fitch_jobs=None,
        max_memory_mb=None,
        **kwargs,
    ):
        if newick_tree is None and species_list is None:
            raise ValueError("Either newick_tree or species_list must be provided.")
        self.species_list = species_list
        self.newick_tree = newick_tree
        self.base_output_dir = base_output_dir
        self.run_id = run_id or "multi_species_run"
        self.output_dir = os.path.join(self.base_output_dir, self.run_id)
        self.aligner_name=aligner_name
        self.aligner_cmd=aligner_cmd
        self.no_cache = no_cache
        self.verbose = verbose

        self.requested_cores = cores
        self.cores = parallel.resolve_cores(cores)
        self.no_parallel = no_parallel
        self.align_jobs = align_jobs
        self.scan_jobs = scan_jobs
        self.fitch_jobs = fitch_jobs
        self.max_memory_mb = max_memory_mb

        kwargs["cores"] = self.cores
        self.params = kwargs

        self.outgroup_name = outgroup
        self.tree = None
        self.terminal_mapping = None
        self.species_dict = {}
        self.reference = None
        self.genomes = {}
        self.alignments = []
        self.pileup_path = None

        os.makedirs(self.output_dir, exist_ok=True)

    @property
    def parallel(self):
        """True when this run should use the parallel stages at all."""
        return (not self.no_parallel) and self.cores > 1

    def _stage_parallel(self, jobs):
        """True when one stage should be parallel (``jobs=1`` switches it off)."""
        return self.parallel and (jobs is None or jobs > 1)

    def run(self):
        log("Starting multi-species mutation extraction pipeline...", self.verbose)
        # Checked up front to avoid long runs that eventually fail without PHYLIP.
        if not check_phylip_available('dnapars'):
            raise RuntimeError(
                "PHYLIP is required for multi-species phylogenetic reconstruction but was not found.\n"
                "Please install PHYLIP via conda: `conda install -c bioconda phylip`"
            )

        if self.newick_tree:
            self.parse_and_annotate_tree()
        else:
            self.parse_and_annotate_list()
        self.download_index_and_fragment()
        self.align_species_to_outgroup()
        self.generate_pileup()
        self._extract_mutations()
        self._reconstruct_phylogeny()
        log("Pipeline completed successfully.", self.verbose)

    def parse_and_annotate_tree(self):
        accession_lookup, default_outgroup = parse_species_accession_from_newick(self.newick_tree)
        if not self.outgroup_name:
            self.outgroup_name = default_outgroup
        self.tree, self.terminal_mapping, self.species_list = annotate_tree_with_indices(self.newick_tree, self.outgroup_name, verbose=self.verbose)

        # Rebuild species_dict in the same order as species_list (outgroup first),
        # so self.genomes and self.alignments follow the same ordering.
        self.species_dict = {name: accession_lookup[name] for name in self.species_list}

        tree_path = os.path.join(self.output_dir, "annotated_tree.nwk")
        save_annotated_tree(self.tree, tree_path)
        with open(os.path.join(self.output_dir, "species_mapping.json"), 'w') as f:
            json.dump(self.terminal_mapping, f, indent=2)
        #with open(os.path.join(self.output_dir, "species_mapping2.json"), 'w') as f:
        #    json.dump(self.species_dict, f, indent=2)

    def parse_and_annotate_list(self):
        if not self.outgroup_name:
            raise ValueError("Outgroup name must be provided when species_list is used.")

        accession_lookup = {key: value for key, value in self.species_list}

        self.species_list, self.terminal_mapping = annotate_list_with_indices(self.species_list, self.outgroup_name, verbose=self.verbose)

        # Rebuild species_dict in the same order as species_list (outgroup first),
        # so self.genomes and self.alignments follow the same ordering.
        self.species_dict = {name: accession_lookup[name] for name in self.species_list}

        with open(os.path.join(self.output_dir, "species_mapping.json"), 'w') as f:
            json.dump(self.terminal_mapping, f, indent=2)
        #with open(os.path.join(self.output_dir, "species_mapping2.json"), 'w') as f:
        #    json.dump(self.species_dict, f, indent=2)


    def download_index_and_fragment(self):
        for species, accession in self.species_dict.items():
            genome = Genome(
                name=species,
                accession=accession,
                output_dir=self.output_dir,
                no_cache=self.no_cache,
                verbose=self.verbose
            )
            genome.download()

            if species == self.outgroup_name:
                genome.index(aligner=self.aligner_name)
                self.reference = genome
            else:
                genome.generate_fragment_fastq(
                    length=self.params.get("fragment_length", 150),
                    offset=self.params.get("fragment_offset", 75),
                    force=self.no_cache
                )
            self.genomes[species] = genome

        with open(os.path.join(self.output_dir, "genome_species_mapping.json"), 'w') as f:
            json.dump(self.species_dict, f, indent=2)

    def align_species_to_outgroup(self):
        ingroup = [(species, genome) for species, genome in self.genomes.items()
                   if species != self.outgroup_name]
        if self._stage_parallel(self.align_jobs) and len(ingroup) > 1:
            return self._align_species_parallel(ingroup)

        for species, genome in self.genomes.items():
            if species == self.outgroup_name:
                continue

            aligner = Aligner(
                species_genome=genome,
                reference_genome=self.reference,
                base_output_dir=self.output_dir,
                aligner_cmd=self.aligner_cmd,
                aligner_name=self.aligner_name,
                cores=self.params.get('cores', None),
                verbose=self.verbose
            )

            if self.params.get("streamed", False):
                aligner.align_streamed(
                    mapq=self.params.get("mapq", 60),
                    low_mapq=self.params.get("low_mapq", 1),
                    max_sort_mem=self.params.get("max_samtools_mem", None),
                    continuity=self.params.get("continuity", True)
                )
            else:
                aligner.align_disk_cached(
                    mapq=self.params.get("mapq", 60),
                    low_mapq=self.params.get("low_mapq", 1),
                    continuity=self.params.get("continuity", True)
                )

            self.alignments.append(aligner)

    def _align_species_parallel(self, ingroup):
        """Align ingroup species concurrently, splitting the CPU budget across species."""
        n_workers, thread_counts = parallel.plan_align_jobs(
            len(ingroup), self.cores, self.align_jobs, self.max_memory_mb)

        # Same construction as the serial stage, with the thread budget split.
        aligners = [
            Aligner(
                species_genome=genome,
                reference_genome=self.reference,
                base_output_dir=self.output_dir,
                aligner_cmd=self.aligner_cmd,
                aligner_name=self.aligner_name,
                cores=n_threads,
                verbose=self.verbose,
            )
            for (_, genome), n_threads in zip(ingroup, thread_counts)
        ]

        streamed = self.params.get("streamed", False)
        align_kwargs = {
            "mapq": self.params.get("mapq", 60),
            "low_mapq": self.params.get("low_mapq", 1),
            "continuity": self.params.get("continuity", True),
        }
        if streamed:
            align_kwargs["max_sort_mem"] = self.params.get("max_samtools_mem", None)

        # `self.genomes` follows `species_list` order, so `aligners` does too
        self.alignments = parallel.run_alignments(
            aligners, streamed, align_kwargs, n_workers, verbose=self.verbose)

    def generate_pileup(self):
        pileup = Pileup(
            outgroup=self.reference,
            aligners=self.alignments,
            base_output_dir=self.output_dir,
            run_id=self.run_id,
            no_cache=self.no_cache,
            verbose=self.verbose
        )

        if not self._stage_parallel(self.scan_jobs):
            self.pileup_path = pileup.generate()
            return

        # The parallel scan pileups each chromosome itself, straight from the
        # indexed BAMs, so the whole-genome pileup is never read. 
        for path in [pileup.ref_fasta] + [a.final_bam for a in self.alignments]:
            pileup._check_file(path)
        self.pileup_path = pileup.pileup_path
        log("Parallel extraction enabled: skipping whole-genome pileup "
            "(per-chromosome pileups are generated during extraction).", self.verbose)


    def _extract_mutations(self):
        # The reference and BAMs are passed only when the scan is parallel.
        parallel_scan_kwargs = {}
        if self._stage_parallel(self.scan_jobs):
            parallel_scan_kwargs = dict(
                ref_fasta=self.reference.fasta_path,
                bams=[a.final_bam for a in self.alignments],
                fai_path=self.reference.fasta_path + ".fai",
            )

        extractor = MultipleSpeciesMutationExtractor(
        pileup_file=self.pileup_path,
        output_dir=self.output_dir,
        n_species=len(self.genomes),
        tree=self.tree,
        species_list=self.species_list,
        mapping=self.terminal_mapping,
        no_cache=False,
        verbose=True,
        cores=self.cores,
        scan_jobs=self.scan_jobs,
        max_memory_mb=self.max_memory_mb,
        fitch_jobs=self.fitch_jobs,
        parallel_fitch=self._stage_parallel(self.fitch_jobs),
        **parallel_scan_kwargs
        )
        extractor.extract()


    def _reconstruct_phylogeny(self):
        run_phylip(
            command='dnapars',
            df_path=os.path.join(self.output_dir, "matching_bases.csv.gz"),
            tree_path=os.path.join(self.output_dir, "annotated_tree.nwk") if self.newick_tree else None,
            output_dir=self.output_dir,
            prefix="multi_species_phylip",
            input_string="5\nY\n",
            mapping=self.terminal_mapping,
            verbose=self.verbose
        )


if __name__ == "__main__":
    species = [
        ("Drosophila_pseudoobscura", "GCF_009870125.1"),
        ("Drosophila_miranda", "GCF_003369915.1")
    ]
    outgroup = ("Drosophila_helvetica", "GCA_963969585.1")

    pipeline = MutationExtractionPipeline(
        species_list=species,
        outgroup=outgroup,
        aligner_name="bwa-mem2",
        base_output_dir="../Output_OO",
        mapq=60, 
        suffix= 'MAPQ60',
        cores=16
    )
    pipeline.run()

    species_list = [("Drosophila_pseudoobscura", "GCF_009870125.1"),
                    ("Drosophila_miranda", "GCF_003369915.1"),
                    ("Drosophila_helvetica", "GCA_963969585.1")]
    run_id = 'drosophila1_run_mutiple_species'
    pipeline = MultiSpeciesMutationPipeline(species_list=species_list,
                                            base_output_dir="../Output_OO",
                                            run_id=run_id,
                                            outgroup='Drosophila_helvetica',
                                            cores=16)
    
    pipeline.run()
    """
    newick_tree = "(((Drosophila_sechellia|GCF_004382195.2,Drosophila_melanogaster|GCF_000001215.4),Drosophila_mauritiana|GCF_004382145.1),Drosophila_santomea|GCF_016746245.2);"
    
    run_id = 'drosophila2_run_mutiple_species'
    pipeline = MultiSpeciesMutationPipeline(newick_tree,
                                            base_output_dir="../Output_OO",
                                            run_id=run_id,
                                            outgroup='Drosophila_santomea')
    """
    """
    species_list = [('Drosophila_sechellia','GCF_004382195.2'),
                    ('Drosophila_melanogaster','GCF_000001215.4'),
                    ('Drosophila_mauritiana', 'GCF_004382145.1'), 
                    ('Drosophila_santomea','GCF_016746245.2')]

    run_id = 'drosophila1_run_mutiple_species'
    pipeline = MultiSpeciesMutationPipeline(species_list=species_list,
                                            base_output_dir="../Output_OO",
                                            run_id=run_id,
                                            outgroup='Drosophila_santomea')
    
    pipeline.run()
    """

