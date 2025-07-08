#!/usr/bin/env python3

from aastk.cli import get_main_parser
from aastk.pasr import *
from aastk.cugo import *
from aastk.casm import *

def main():
    parser = get_main_parser()
    args = parser.parse_args()

    # If no subcommand is specified, print help
    if not args.subparser_name:
        parser.print_help()
        return

    try:
        ### PARSER FOR PASR FUNCTIONALITIES AND WORKFLOW ###
        if args.subparser_name == 'build':
            db_path = build_protein_db(
                db_dir=args.db,
                protein_name=args.protein_name,
                seed_fasta=args.seed,
                threads=args.threads,
                force=args.force
            )

        elif args.subparser_name == 'search':
            search_results = search_protein_db(
                db_path=args.db,
                query_path=args.query,
                protein_name=args.protein_name,
                output_dir=args.output,
                threads=args.threads,
                sensitivity=args.sensitivity,
                block=args.block,
                chunk=args.chunk,
                force=args.force
            )

        elif args.subparser_name == 'extract':
            extract_matching_sequences(
                protein_name=args.protein_name,
                blast_tab=args.tabular,
                query_path=args.query,
                output_dir=args.output,
                key_column=args.key_column,
                force=args.force
            )

        elif args.subparser_name == 'calculate':
            calculate_max_scores(
                protein_name=args.protein_name,
                extracted=args.extracted,
                matrix=args.matrix,
                output_dir=args.output
            )

        elif args.subparser_name == 'bsr':
            blast_score_ratio(
                protein_name=args.protein_name,
                blast_tab=args.tabular,
                max_scores_path=args.max_scores,
                output_dir=args.output,
                key_column=args.key_column,
                column_info_path=args.column_info_path,
                score_column=args.score_column,
                force=args.force
            )

        elif args.subparser_name == 'pasr_plot':
            plot_bsr(
                protein_name=args.protein_name,
                bsr_file=args.bsr,
                output_dir=args.output,
                yaml_path=args.yaml,
                force=args.force,
                update=args.update
            )

        elif args.subparser_name == 'metadata':
            metadata(
                selfmin=args.selfmin,
                selfmax=args.selfmax,
                dbmin=args.dbmin,
                bsr=args.bsr_cutoff,
                output_dir=args.output,
                protein_name=args.protein_name,
                force=args.force
            )

        elif args.subparser_name == 'select':
            select(
                yaml_path=args.yaml,
                matched_fasta=args.matched,
                bsr_table=args.bsr,
                output_dir=args.ouput,
                selfmin=args.selfmin,
                selfmax=args.selfmax,
                dbmin=args.dbmin,
                bsr=args.bsr_cutoff,
                protein_name=args.protein_name,
                force=args.force,
                create_yaml=args.create_yaml,
                params=args.params
            )

        elif args.subparser_name == 'pasr':
            pasr(
                protein_name=args.protein_name,
                seed_fasta=args.seed,
                query_fasta=args.query,
                matrix_name=args.matrix,
                threads=args.threads,
                output_dir=args.output,
                block=args.block,
                chunk=args.chunk,
                sensitivity=args.sensitivity,
                update=args.update,
                yaml_path=args.yaml,
                force=args.force
            )

        ### PARSER FOR CUGO FUNCTIONALITIES AND WORKFLOW ###
        elif args.subparser_name == 'parse':
            parse(
                tar_gz_path=args.gff_path,
                tmhmm_tar_path=args.tmhmm_dir,
                output_dir=args.output,
                globdb_version=args.globdb_version,
                force=args.force
            )


        elif args.subparser_name == 'context':
            context(
                protein_ids=args.protein_ids,
                cugo_path=args.cugo_path,
                cugo_range=args.cugo_range,
                output_dir=args.output,
                protein_name=args.protein_name,
                force=args.force,
                fasta_path=args.fasta
            )

        elif args.subparser_name == 'cugo_plot':
            cugo_plot(
                cugo_path=args.cugo_path,
                flank_lower=args.flank_lower,
                flank_upper=args.flank_upper,
                top_n=args.top_n,
                cugo=args.cugo,
                size=args.size,
                all_plots=args.all,
                bin_width=args.bin_width,
                y_range=args.y_range
            )

        ### PARSER FOR ASM_CLUST FUNCTIONALITIES AND WORKFLOW ###
        elif args.subparser_name == 'casm':
            casm(
                fasta=args.fasta,
                output=args.output,
                subset=args.subset,
                subset_size=args.subset_size,
                threads=args.threads,
                perplexity=args.perplexity,
                iterations=args.iterations,
                exaggeration=args.exaggeration,
                metadata_protein=args.metadata_protein,
                metadata_genome=args.metadata_genome,
                force=args.force,
            )

        elif args.subparser_name == 'casm_plot':
            casm_plot(
                tsv_file=args.tsv,
                output=args.output
            )
    except Exception as e:
        logger.error(f"Error executing command: {e}")
        return 1
    return 0

if __name__ == "__main__":
    exit(main())