# species_mutation_extraction/mutextractor/pipeline.py

import json
import os
import time 
import gc

import pandas as pd
from .cleanup_manager import PipelineCleaner
from .genome_manager import Genome
from .alignment_manager import Aligner
from .multiple_species_mutation_extractor_manager import MultipleSpeciesMutationExtractor
from .mutation_extractor_manager import FiveMerExtractor, MutationExtractor, MutationNormalizer, TripletExtractor
from .pileup_manager import Pileup
from .plot_utils import CoveragePlotter, MutationDensityPlotter, MutationSpectraPlotter
from .utils import get_top_n_chromosomes, log
import psutil
import pysam

class MutationExtractionPipeline:
    def __init__(self, 
                 species_list,
                 outgroup,
                 aligner_name="bwa", 
                 aligner_cmd=None,
                 base_output_dir="../Output", 
                 no_cache = False,
                 verbose = True, 
                 run_id = None,
                 **kwargs):
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
        self.verbose = verbose
        self.no_cache = no_cache

    
    def run(self):
        log("Starting mutation extraction pipeline...", self.verbose)
        timings = {}
        memory_log = {}
        process = psutil.Process(os.getpid())
        start_pipeline = time.time()

        def get_memory():
            return round(process.memory_info().rss / (1024 ** 2), 2)  # In MB

        def timed_stage(stage_name, func):
            log(f"--- Starting: {stage_name} ---", self.verbose)
            mem_before = get_memory()
            start = time.time()
            func()
            gc.collect()  # Clean up memory after each stage
            end = time.time()
            mem_after = get_memory()

            timings[stage_name] = round(end - start, 2)
            memory_log[stage_name] = {"start_MB": mem_before, "end_MB": mem_after}
            log(f"{stage_name} completed in {timings[stage_name]} seconds", self.verbose)
            log(f"Memory usage: {mem_before} → {mem_after} MB", self.verbose)

        timed_stage("Download and Fragment Genomes", self.download_index_and_fragment_genomes)
        timed_stage("Align Species", self.align_species)
        timed_stage("Generate Pileup", self.generate_pileup)
        timed_stage("Extract Mutations and Triplets", self.extract_mutations_and_triplets)
        timed_stage("Extract Intervals", self.extract_intervals)
        timed_stage("Run Plots", self.run_plots)
        timed_stage("Cleanup files", self.cleanup)

        total_runtime = round(time.time() - start_pipeline, 2)
        timings["Total Runtime"] = total_runtime
        memory_log["Total Runtime"] = {"final_MB": get_memory()}

        timing_path = os.path.join(self.output_dir, "pipeline_timings.json")
        with open(timing_path, "w") as f:
            json.dump({"timings": timings, "memory": memory_log}, f, indent=2)

        log(f"Timing and memory info saved to: {timing_path}", self.verbose)
        log("Pipeline completed successfully.", self.verbose)


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
            genome.generate_fragment_fastq(length=self.params.get("fragment_length", 150), offset=self.params.get("fragment_offset", 75), force=self.no_cache)
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
        
        self.pileup_path = pileup_generator.generate()


    def extract_mutations_and_triplets(self):
        # log("Extracting 3mer mutations and triplets from pileup...", self.verbose)
        mutation_extractor = MutationExtractor(reference=self.reference.name,
                              taxon1=self.genomes[0].name,
                              taxon2=self.genomes[1].name,
                              pileup_file=self.pileup_path,
                              mutation_output_dir=os.path.join(self.output_dir, 'Mutations'),
                              triplet_output_dir=os.path.join(self.output_dir, 'Triplets'),
                              no_full_mutations=False,
                              no_cache=False,
                              verbose=self.verbose)
        mutation_extractor.extract()

        fivemer_extractor = FiveMerExtractor(reference=self.reference.name,
                              taxon1=self.genomes[0].name,
                              taxon2=self.genomes[1].name,
                              pileup_file=self.pileup_path,
                              output_dir=os.path.join(self.output_dir, 'Mutations'),
                              no_cache=False,
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
            verbose=True,
            divergence_time= self.params.get("divergence_time", None),
        )
        normalizer.normalize()

    def _extract_bam_intervals(self, input_bam, output_dir, sorted=False, merge=False, no_cache=False):
            os.makedirs(output_dir, exist_ok=True)

            base_name = os.path.basename(input_bam).rsplit(".", 1)[0]
            output_file = os.path.join(output_dir, f"{base_name}_intervals.tsv.gz")

            if os.path.exists(output_file) and not no_cache:
                log(f"Intervals already exist: {output_file}", self.verbose)
                return output_file

            bamfile = pysam.AlignmentFile(input_bam, "rb")

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

            intervals = (
                extract_intervals_sorted(bamfile)
                if sorted
                else merge_intervals(extract_raw_intervals(bamfile)) if merge
                else extract_raw_intervals(bamfile)
            )

            df = pd.DataFrame(intervals, columns=["chromosome", "start", "end"])
            df.to_csv(output_file, sep='\t', index=False, compression="gzip")

            log(f"Intervals written to: {output_file}", self.verbose)
            return output_file
    
    def extract_intervals(self):
        for bam in self.alignments:
            self._extract_bam_intervals(bam.final_bam, os.path.join(self.output_dir, 'Intervals'))    

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
        cleaner = PipelineCleaner(self.genomes + [self.reference], self.alignments, self.pileup, base_dir=self.output_dir, verbose=True)
        cleaner.run(bams=True, pileup=True, genomes=True)


from .multiple_species_utils import (
    annotate_list_with_indices,
    parse_species_accession_from_newick,
    annotate_tree_with_indices,
    save_annotated_tree,
)
from .run_phylip import run_phylip, check_phylip_available


class MultiSpeciesMutationPipeline:
    def __init__(
        self,
        newick_tree = None,
        species_list = None,
        base_output_dir="../Output",
        run_id=None,
        outgroup=None,
        aligner_name="bwa", 
        aligner_cmd=None,
        no_cache=False,
        verbose=True,
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

    def run(self):
        log("Starting multi-species mutation extraction pipeline...", self.verbose)
        if self.newick_tree:
            self.parse_and_annotate_tree()
        else:
            self.parse_and_annotate_list()
        self.download_index_and_fragment()
        self.align_species_to_outgroup()
        self.generate_pileup()
        self._extract_mutations()
        
        # Validate PHYLIP is available before phylogenetic reconstruction
        if not check_phylip_available('dnapars'):
            raise RuntimeError(
                "PHYLIP is required for multi-species phylogenetic reconstruction but was not found.\n"
                "Please install PHYLIP via conda: `conda install -c bioconda phylip`"
            )
        
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

    def generate_pileup(self):
        pileup = Pileup(
            outgroup=self.reference,
            aligners=self.alignments,
            base_output_dir=self.output_dir,
            run_id=self.run_id,
            no_cache=self.no_cache,
            verbose=self.verbose
        )
        self.pileup_path = pileup.generate()


    def _extract_mutations(self):
        extractor = MultipleSpeciesMutationExtractor(
        pileup_file=self.pileup_path,
        output_dir=self.output_dir,
        n_species=len(self.genomes),
        tree=self.tree,
        species_list=self.species_list,
        mapping=self.terminal_mapping,
        no_cache=False,
        verbose=True
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
        aligner="bwa",
        base_output_dir="../Output_OO",
        mapq=60, 
        suffix= 'MAPQ60'
    )
    pipeline.run()

    # newick_tree = "(((Drosophila_sechellia|GCF_004382195.2,Drosophila_melanogaster|GCF_000001215.4),Drosophila_mauritiana|GCF_004382145.1),Drosophila_santomea|GCF_016746245.2);"
    """
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

