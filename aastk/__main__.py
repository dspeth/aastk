#!/usr/bin/env python3

import logging
from .log import logger_setup
from aastk.cli import get_main_parser
from aastk.pasr import *
from aastk.cugo import *
from aastk.casm import *

def main():
    parser = get_main_parser()
    args = parser.parse_args()

    logger = logger_setup(silent=args.silent)

    # If no subcommand is specified, print help
    if not args.subparser_name:
        parser.print_help()
        return

    try:
        ### PARSER FOR PASR FUNCTIONALITIES AND WORKFLOW ###
        if args.subparser_name == 'build':
            build(
                db_dir=args.db,
                seed_fasta=args.seed,
                threads=args.threads,
                force=args.force
            )

        elif args.subparser_name == 'search':
            search(
                db_path=args.db,
                query_path=args.query,
                output_dir=args.output,
                threads=args.threads,
                sensitivity=args.sensitivity,
                block=args.block,
                chunk=args.chunk,
                force=args.force
            )

        elif args.subparser_name == 'extract':
            extract(
                blast_tab=args.tabular,
                query_path=args.query,
                output_dir=args.output,
                key_column=args.key_column,
                force=args.force
            )

        elif args.subparser_name == 'calculate':
            calculate(
                extracted=args.extracted,
                matrix=args.matrix,
                output_dir=args.output
            )

        elif args.subparser_name == 'bsr':
            bsr(
                blast_tab=args.tabular,
                max_scores_path=args.max_scores,
                output_dir=args.output,
                key_column=args.key_column,
                column_info_path=args.column_info_path,
                score_column=args.score_column,
                force=args.force
            )

        elif args.subparser_name == 'pasr_plot':
            pasr_plot(
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
                force=args.force
            )

        elif args.subparser_name == 'select':
            select(
                yaml_path=args.yaml,
                matched_fasta=args.matched,
                bsr_table=args.bsr,
                output_dir=args.output,
                selfmin=args.selfmin,
                selfmax=args.selfmax,
                dbmin=args.dbmin,
                bsr=args.bsr_cutoff,
                force=args.force,
                create_yaml=args.create_yaml,
                params=args.params
            )

        elif args.subparser_name == 'pasr':
            pasr(
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
                cog_gff_tar_path=args.cog_gff,
                kegg_gff_tar_path=args.kegg_gff,
                pfam_gff_tar_path=args.pfam_gff,
                tmhmm_tar_path=args.tmhmm_dir,
                taxonomy_path=args.taxonomy_path,
                output_dir=args.output,
                globdb_version=args.globdb_version,
                force=args.force
            )


        elif args.subparser_name == 'context':
            context(
                fasta=args.fasta,
                id_list=args.id_list,
                cugo_path=args.cugo_path,
                output_dir=args.output,
                threads=args.threads,
                force=args.force,
            )

        elif args.subparser_name == 'cugo_plot':
            cugo_plot(
                context_path=args.context_path,
                flank_lower=args.flank_lower,
                flank_upper=args.flank_upper,
                top_n=args.top_n,
                output=args.output,
                cugo=args.cugo,
                size=args.size,
                all_plots=args.all,
                bin_width=args.bin_width,
                y_range=args.y_range,
                tmh_y_range=args.tmh_y_range,
                force=args.force
            )

        elif args.subparser_name == 'cugo':
            cugo(
                cugo_path=args.cugo_path,
                fasta=args.fasta,
                id_list=args.id_list,
                output_dir=args.output,
                flank_lower=args.flank_lower,
                flank_upper=args.flank_upper,
                top_n=args.top_n,
                threads=args.threads,
                force=args.force,
                bin_width=args.bin_width,
                y_range=args.y_range,
                tmh_y_range=args.tmh_y_range
            )

        elif args.subparser_name == 'retrieve':
            retrieve(
                context_path=args.context_path,
                position=args.position,
                output=args.output
            )

        ### PARSER FOR CASM FUNCTIONALITIES AND WORKFLOW ###
        elif args.subparser_name == 'matrix':
            matrix(
                fasta=args.fasta,
                output=args.output,
                subset=args.subset,
                subset_size=args.subset_size,
                threads=args.threads,
                force=args.force
            )

        elif args.subparser_name == 'cluster':
            cluster(
                matrix_path=args.matrix_path,
                matrix_metadata_path=args.metadata_matrix,
                output=args.output,
                perplexity=args.perplexity,
                iterations=args.iterations,
                exaggeration=args.exagerration,
                threads=args.threads,
                metadata_protein=args.metadata_protein,
                metadata_genome=args.metadata_genome,
                force=args.force,
            )


        elif args.subparser_name == 'casm_plot':
            casm_plot(
                early_clust_path=args.early_clust,
                full_clust_path=args.full_clust,
                output=args.output,
                show_cluster_numbers=args.show
            )

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
                show_cluster_numbers=args.show
            )

        elif args.subparser_name == 'pick':
            pick(
                final_embedding_file=args.full_clust,
                fasta=args.fasta,
                no_cluster=args.no_cluster,
                output=args.output,
                force=args.force
            )


    except Exception as e:
        logger.error(f"Error executing command: {e}")
        return 1
    return 0

if __name__ == "__main__":
    exit(main())
