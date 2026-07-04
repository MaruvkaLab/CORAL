import argparse
import json
import os
import sys
from .pipeline import MutationExtractionPipeline, MultiSpeciesMutationPipeline
from .run_phylip import run_phylip

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
    single.add_argument("--aligner-name", default="bwa")
    single.add_argument("--aligner-cmd", default=None)
    single.add_argument("--streamed", action="store_true")
    single.add_argument("--mapq", type=int, default=60)
    single.add_argument("--low-mapq", type=int, default=1)
    continuity_group = single.add_mutually_exclusive_group()
    continuity_group.add_argument("--continuity", action="store_true", default=True, help="Enable continuity mode (default)")
    continuity_group.add_argument("--no-continuity", dest="continuity", action="store_false", help="Disable continuity mode")
    single.add_argument("--cores", type=int, default=None)
    single.add_argument("--divergence-time", type=int, default=None)
    single.add_argument("--five-mer", dest="five_mer", action="store_true",
                        help="Also extract 5-mer mutation contexts (adds a full pileup pass; off by default)")

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
    multi.add_argument("--aligner-name", default="bwa")
    multi.add_argument("--aligner-cmd", default=None)
    multi.add_argument("--streamed", action="store_true")
    continuity_group_multi = multi.add_mutually_exclusive_group()
    continuity_group_multi.add_argument("--continuity", action="store_true", default=True, help="Enable continuity mode (default)")
    continuity_group_multi.add_argument("--no-continuity", dest="continuity", action="store_false", help="Disable continuity mode")
    multi.add_argument("--run-id", default=None)
    multi.add_argument("--mapq", type=int, default=60)
    multi.add_argument("--low-mapq", type=int, default=1)
    multi.add_argument("--cores", type=int, default=None)

    # === Run PHYLIP ===
    phylip = subparsers.add_parser("run_phylip", help="Run PHYLIP on mutation matrix")
    phylip.add_argument("--df", required=True, help="Path to matching_bases.csv.gz")
    phylip.add_argument("--tree", required=True, help="Path to annotated_tree.nwk")
    phylip.add_argument("--mapping", required=True, help="Path to species_mapping.json")
    phylip.add_argument("--phylip-command", dest="phylip_command", default="dnapars", help="PHYLIP command: dnapars, dnapenny, etc.")
    phylip.add_argument("--prefix", default="phylip_run")
    phylip.add_argument("--input-string", default="Y\n")
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
            )
            pipeline.run()

        elif args.subcmd == "run_multi":
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
            )
            pipeline.run()

        elif args.subcmd == "run_phylip":
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
                verbose=args.verbose
            )

        else:
            parser.print_help()
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()

