import argparse
import json
import os
import sys
from .pipeline import MutationExtractionPipeline, MultiSpeciesMutationPipeline

def main():
    parser = argparse.ArgumentParser(description="Species Mutation Extraction CLI")
    subparsers = parser.add_subparsers(dest="subcmd")

    # === Single-Pipeline ===
    single = subparsers.add_parser("run_single", help="Run 3-species pipeline (outgroup + 2)")
    single.add_argument("--outgroup", nargs=2, metavar=("NAME", "ACCESSION"), required=True)
    single.add_argument("--species", nargs=4, metavar=("NAME1", "ACC1", "NAME2", "ACC2"), required=True)
    single.add_argument("--output", required=True)
    single.add_argument("--no-cache", action="store_true")
    verbose_group = single.add_mutually_exclusive_group()
    verbose_group.add_argument("--verbose", dest="verbose", action="store_true", help="Enable verbose logging (default: enabled)")
    verbose_group.add_argument("--quiet", dest="verbose", action="store_false", help="Disable verbose logging")
    single.set_defaults(verbose=True)
    single.add_argument("--suffix", default=None)
    single.add_argument("--aligner-name", default="bwa-mem2")
    single.add_argument("--aligner-cmd", default=None)
    single.add_argument("--streamed", action="store_true")
    single.add_argument("--mapq", type=int, default=60)
    single.add_argument("--low-mapq", type=int, default=1)
    continuity_group = single.add_mutually_exclusive_group()
    continuity_group.add_argument("--continuity", action="store_true", default=True, help="Enable continuity mode (default)")
    continuity_group.add_argument("--no-continuity", dest="continuity", action="store_false", help="Disable continuity mode")
    single.add_argument("--cores", type=int, default=None)
    single.add_argument("--divergence-time", type=int, default=None)
    single.add_argument("--no-plots", dest="plots", action="store_false", default=True,
                        help="Skip the plotting stage (saves time and a lot of small PNG files)")
    single.add_argument("--five-mer", dest="five_mer", action="store_true",
                        help="Also extract 5-mer mutation contexts (adds a full pileup pass; off by default)")
    single.add_argument("--repeat-mask", dest="repeat_mask", nargs="?", const="windowmasker",
                        default=None, choices=["windowmasker", "repeatmasker"], metavar="TOOL",
                        help="Mask repeats: skip pseudo-reads from species repeats and calls at "
                             "reference repeat positions. Off unless given; TOOL is windowmasker "
                             "(default, library-free) or repeatmasker")
    single.add_argument("--repeat-species", dest="repeat_species", default=None, metavar="CLADE",
                        help="RepeatMasker library clade (e.g. 'drosophila'); only with --repeat-mask repeatmasker")

    multi = subparsers.add_parser("run_multi", help="Run multi-species pipeline from Newick")
    multi.add_argument("--newick-tree", default=None)
    multi.add_argument("--species-list", type=json.loads, default=None, help='List of species as JSON, e.g. \'[["Homo_sapiens", "GCF_..."], ...]\'')
    multi.add_argument("--output", required=True)
    multi.add_argument("--no-cache", action="store_true")
    verbose_group_multi = multi.add_mutually_exclusive_group()
    verbose_group_multi.add_argument("--verbose", dest="verbose", action="store_true", help="Enable verbose logging (default: enabled)")
    verbose_group_multi.add_argument("--quiet", dest="verbose", action="store_false", help="Disable verbose logging")
    multi.set_defaults(verbose=True)
    multi.add_argument("--outgroup", default=None)
    multi.add_argument("--aligner-name", default="bwa-mem2")
    multi.add_argument("--aligner-cmd", default=None)
    multi.add_argument("--streamed", action="store_true")
    continuity_group_multi = multi.add_mutually_exclusive_group()
    continuity_group_multi.add_argument("--continuity", action="store_true", default=True, help="Enable continuity mode (default)")
    continuity_group_multi.add_argument("--no-continuity", dest="continuity", action="store_false", help="Disable continuity mode")
    multi.add_argument("--run-id", default=None)
    multi.add_argument("--mapq", type=int, default=60)
    multi.add_argument("--low-mapq", type=int, default=1)
    multi.add_argument("--cores", type=int, default=None,
                       help="Number of cores to use. Default: use all available cores. "
                       "--cores 1 runs everything serially.")

    # --- parallelism (run_multi is parallel by default; these switch it off) ---
    multi.add_argument("--no-parallel", action="store_true",
                       help="Run the mutation extraction serially: whole-genome pileup, "
                            "serial scan, recursive per-row Fitch.")
    multi.add_argument("--align-jobs", type=int, default=None,
                       help="Species aligned concurrently; each gets cores/jobs threads. ")
    multi.add_argument("--scan-jobs", type=int, default=None,
                       help="Chromosomes scanned concurrently. 1 uses the whole-genome pileup.")
    multi.add_argument("--fitch-jobs", type=int, default=None,
                       help="Workers for the Fitch pass. 1 runs the recursive per-row "
                            "implementation.")
    multi.add_argument("--max-memory-mb", type=int, default=None,
                       help="Memory budget for the whole run. Caps the worker count of "
                            "each stage.")
    multi.add_argument("--skip-no-tree", action="store_true",
                       help="With --newick-tree, skip the PHYLIP search without a starting "
                            "tree. The spectra use the given tree either way; only the "
                            "parsimony score comparison is lost.")
    multi.add_argument("--no-node-sequences", action="store_true",
                       help="With --newick-tree, run PHYLIP on the given tree without menu "
                            "option 5 (print sequences at all nodes). CORAL never reads that "
                            "table, and on large runs it is what runs out of memory; the tree, "
                            "the edge table and the score are unchanged.")
    multi.add_argument("--given-tree-sites", type=int, default=None, metavar="N",
                       help="With --newick-tree, run PHYLIP on the given tree with a seeded "
                            "subsample of N of its sites. CORAL reads only the tree and the "
                            "edge table from that run, which come from the given tree; the "
                            "score changes and is no longer compared. Fitch still uses every site.")

    # --- split jobs: one HTCondor job per species for the alignment (see split_jobs.py) ---
    split = multi.add_argument_group(
        "split jobs (HTCondor)",
        "Off unless --split-jobs is given. Submits the run as a DAG: one job downloads and "
        "indexes the outgroup, one job per ingroup species downloads, fragments and aligns "
        "it, and a final job runs the rest of run_multi on the cached alignments. The output "
        "directory must be on a filesystem the execute nodes share.")
    split.add_argument("--split-jobs", action="store_true",
                       help="Submit the run as an HTCondor DAG instead of running it here.")
    split.add_argument("--job-cores", type=int, default=4,
                       help="request_cpus of each alignment job; the aligner and sort get "
                            "one thread less (default: 4).")
    split.add_argument("--job-memory-mb", type=int, default=8000,
                       help="request_memory of each alignment job and of the outgroup job "
                            "(default: 8000).")
    split.add_argument("--job-extra", action="append", default=None, metavar="LINE",
                       help="Extra submit-description line for the alignment and outgroup "
                            "jobs, e.g. 'request_disk = 20G'. Repeatable; written last, so "
                            "it overrides the defaults.")
    split.add_argument("--main-job-cores", type=int, default=None,
                       help="request_cpus of the final job (default: 2 x --cores).")
    split.add_argument("--main-job-memory-mb", type=int, default=None,
                       help="request_memory of the final job (default: --max-memory-mb, "
                            "else 16384).")
    split.add_argument("--main-job-extra", action="append", default=None, metavar="LINE",
                       help="Extra submit-description line for the final job. Repeatable.")
    split.add_argument("--job-max-concurrent", type=int, default=20,
                       help="Alignment jobs DAGMan runs at once; 0 for no limit (default: 20).")
    split.add_argument("--job-log-dir", default=None,
                       help="Directory for the jobs' .log/.out/.err "
                            "(default: <output>/<run-id>/split_jobs/logs).")
    split.add_argument("--no-submit", action="store_true",
                       help="Write the DAG and submit files but do not submit them.")
    # Set by the generated jobs themselves, not by users.
    split.add_argument("--split-stage", choices=["prepare", "align", "main"], default=None,
                       help=argparse.SUPPRESS)
    split.add_argument("--split-index", type=int, default=None, help=argparse.SUPPRESS)

    # === Run PHYLIP ===
    phylip = subparsers.add_parser("run_phylip", help="Run PHYLIP on mutation matrix")
    phylip.add_argument("--df", required=True, help="Path to matching_bases.csv.gz")
    phylip.add_argument("--tree", required=True, help="Path to annotated_tree.nwk")
    phylip.add_argument("--mapping", required=True, help="Path to species_mapping.json")
    phylip.add_argument("--phylip-command", dest="phylip_command", default="dnapars", help="PHYLIP command: dnapars, dnapenny, etc.")
    phylip.add_argument("--prefix", default="phylip_run")
    phylip.add_argument("--input-string", default="Y\n")
    phylip.add_argument("--no-cache", action="store_true")
    verbose_group_phylip = phylip.add_mutually_exclusive_group()
    verbose_group_phylip.add_argument("--verbose", dest="verbose", action="store_true", help="Enable verbose logging (default: enabled)")
    verbose_group_phylip.add_argument("--quiet", dest="verbose", action="store_false", help="Disable verbose logging")
    phylip.set_defaults(verbose=True)

    args = parser.parse_args()

    try:
        if args.subcmd == "run_single":
            pipeline = MutationExtractionPipeline(
                species_list=[(args.species[0], args.species[1]), (args.species[2], args.species[3])],
                outgroup=(args.outgroup[0], args.outgroup[1]),
                base_output_dir=args.output,
                no_cache=args.no_cache,
                verbose=args.verbose,
                suffix=args.suffix,
                aligner_name=args.aligner_name,
                aligner_cmd=args.aligner_cmd,
                streamed=args.streamed,
                mapq=args.mapq,
                low_mapq=args.low_mapq,
                cores=args.cores,
                continuity=args.continuity,
                divergence_time=args.divergence_time,
                five_mer=args.five_mer,
                plots=args.plots,
                repeat_mask=bool(args.repeat_mask),
                repeat_masker=args.repeat_mask or "windowmasker",
                repeat_species=args.repeat_species,
            )
            pipeline.run()

        elif args.subcmd == "run_multi":
            if args.skip_no_tree and not args.newick_tree:
                parser.error("--skip-no-tree needs --newick-tree")
            if args.no_node_sequences and not args.newick_tree:
                parser.error("--no-node-sequences needs --newick-tree")
            if args.given_tree_sites is not None and not args.newick_tree:
                parser.error("--given-tree-sites needs --newick-tree")
            if args.given_tree_sites is not None and args.given_tree_sites < 1:
                parser.error("--given-tree-sites must be a positive number of sites")
            pipeline = MultiSpeciesMutationPipeline(
                newick_tree=args.newick_tree,
                species_list=args.species_list,
                base_output_dir=args.output,
                outgroup=args.outgroup,
                no_cache=args.no_cache,
                verbose=args.verbose,
                aligner_name=args.aligner_name,
                aligner_cmd=args.aligner_cmd,
                run_id=args.run_id,
                streamed=args.streamed,
                mapq=args.mapq,
                low_mapq=args.low_mapq,
                cores=args.cores,
                continuity=args.continuity,
                no_parallel=args.no_parallel,
                align_jobs=args.align_jobs,
                scan_jobs=args.scan_jobs,
                fitch_jobs=args.fitch_jobs,
                max_memory_mb=args.max_memory_mb,
                skip_no_tree=args.skip_no_tree,
                node_sequences=not args.no_node_sequences,
                given_tree_sites=args.given_tree_sites,
            )
            if args.split_jobs or args.split_stage:
                from .split_jobs import dispatch
                dispatch(pipeline, args, sys.argv[1:])
            else:
                pipeline.run()

        elif args.subcmd == "run_phylip":
            from .run_phylip import run_phylip
            with open(args.mapping) as f:
                mapping = json.load(f)
            run_phylip(
                command=args.phylip_command,
                df_path=args.df,
                tree_path=args.tree,
                output_dir=os.path.dirname(args.df),
                prefix=args.prefix,
                input_string=args.input_string,
                mapping=mapping,
                verbose=args.verbose,
                no_cache=args.no_cache
            )

        else:
            parser.print_help()
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()

