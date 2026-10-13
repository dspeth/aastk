#!/usr/bin/env python3

import logging
from aastk.log import logger_setup
from aastk.cli import get_main_parser
from aastk.pasr import *
from aastk.cugo import *
from aastk.casm import *
from aastk.db.build import database
from aastk.db.check import database_check
from aastk.db.meta import meta, list_metadata, parse_query
from aastk.db.export import export_fasta, export_by_metadata
from aastk.version import __version__, __copyright__, __author__
import sys

def print_help():
    print('''\

            ...::: AASTK v%s :::...

  Usage:
    aastk <command> <arguments>    to run the tools or commands
    aastk <command> -h             for command specific help
    aastk --silent <command>       to suppress all console output except errors

  Main Tools:
    pasr           Generate or update a comprehensive dataset of homologous proteins
                      Runs the subcommands: build, search, get_hit_seqs, max_score, bsr, and pasr_plot
    casm           Cluster and visualize a complex dataset of proteins (eg. a superfamily)
                      Runs the subcommands: matrix, cluster, and casm_plot
    cugo           Retrieve, calculate, and visualise consensus genomic context of protein data sets
                      Runs the subcommands: context, cugo_plot
    meta           Retrieve protein metadata from AASTK SQLite database

  Helper tools:
    pasr_select    Select target sequences from pasr run based on bsr and score cutoffs
    casm_select    Select cluster(s) from casm analysis, and retrieve sequences
    cugo_select    Select genomic position from cugo analysis and retrieve sequences
    filter         Filter sequence dataset to remove non-homologous sequences
    list_metadata  List available metadata

  Plotting tools:
    pasr_plot      Scatterplot with max score on the x-axis and score against the seed db on y-axis
    casm_plot      Scatterplot of sequences with tSNE coordinates as axes
    cugo_plot      Consensus genomic context plot of annotation, length, and transmembrane segments

  Database utilities
    database       Construct an SQL database from amino acid fasta sequences and associated metadata
    database_check Check whether the SQL database is internally consistent
    export_fasta   Export a fasta file from the AASTK SQL database
    export         Export sequences from the AASTK SQL database matching a metadata query
    ''' % __version__)


def main():
    parser = get_main_parser()

    if len(sys.argv) == 1:
        print_help()
        sys.exit(0)
    elif sys.argv[1] in {'-v', '--v', '-version', '--version'}:
        print(f"AASTK: version {__version__} {__copyright__} {__author__}")
        sys.exit(0)
    elif sys.argv[1] in {'-h', '--h', '-help', '--help'}:
        print_help()
        sys.exit(0)
    else:
        args = parser.parse_args()

    output_dir = getattr(args, 'output', getattr(args, 'output_dir', None))
    logger = logger_setup(silent=args.silent, output_dir=output_dir)

    # If no subcommand is specified, print help
    if not args.subparser_name:
        parser.print_help()
        return

    try:
        ### PARSER FOR PASR FUNCTIONALITIES AND WORKFLOW ###
        if args.subparser_name == 'pasr_plot':
            pasr_plot(
                bsr_file=args.bsr,
                output_dir=args.output,
                yaml_path=args.yaml,
                svg=args.svg,
                force=args.force,
                update=args.update,
                density=args.density
            )

        elif args.subparser_name == 'pasr_select':
            pasr_select(
                yaml_path=args.yaml,
                matched_fasta=args.matched,
                bsr_table=args.bsr,
                output_dir=args.output,
                max_score_min=args.max_score_min,
                max_score_max=args.max_score_max,
                dbmin=args.dbmin,
                bsr=args.bsr_cutoff,
                threads=args.threads,
                force=args.force,
                params=args.params,
                filter_seqs=args.filter_seqs,
                density=args.density
            )

        elif args.subparser_name == 'pasr':
            pasr(
                seed_fasta=args.seed,
                query_fasta=args.query,
                matrix=args.matrix,
                threads=args.threads,
                output_dir=args.output,
                db_path=args.db_path,
                key_column=args.key_column,
                block=args.block,
                chunk=args.chunk,
                sensitivity=args.sensitivity,
                update=args.update,
                yaml_path=args.yaml,
                sql=args.sql,
                keep=args.keep,
                svg=args.svg,
                density=args.density,
                force=args.force
            )

        ### PARSER FOR CUGO FUNCTIONALITIES AND WORKFLOW ###
        elif args.subparser_name == 'database':
            database(
                chunk_dir=args.chunk_dir,
                db_version=args.db_version,
                output=args.output,
                tmhmm_tar_path=args.tmhmm_tar_path,
                taxonomy_path=args.taxonomy_path,
                culture_collection_path=args.culture_collection_path,
                high_level_environment_path=args.high_level_environment_path,
                low_level_environment_path=args.low_level_environment_path,
                anvio=args.anvio,
                tmhmm=args.tmhmm,
                taxonomy=args.taxonomy,
                culture_collection=args.culture_collection,
                environmental_data=args.environmental_data,
                all_sources=args.all_sources,
                threads=args.threads,
                force=args.force,
            )

        elif args.subparser_name == 'database_check':
            database_check(
                db_path=args.db_path,
                output=args.output,
                force=args.force
            )

        elif args.subparser_name == 'meta':
            meta(
                db_path=args.db_path,
                fasta=args.fasta,
                id_list=args.id_list,
                output=args.output,
                threads=args.threads,
                include_annotation=args.include_annotation,
                include_taxonomy=args.include_taxonomy,
                include_culture_collection=args.include_culture_collection,
                include_high_level_environment=args.include_high_level_environment,
                include_low_level_environment=args.include_low_level_environment,
                include_high_level_env_category=args.include_high_level_env_category,
                include_low_level_env_category=args.include_low_level_env_category,
                all_metadata=args.all_metadata,
                force=args.force
            )

        elif args.subparser_name == 'cugo_plot':
            cugo_plot(
                context_path=args.context_path,
                flank_lower=args.flank_lower,
                flank_upper=args.flank_upper,
                annotation=args.annotation,
                top_n=args.top_n,
                output=args.output,
                cugo=args.cugo,
                size=args.size,
                tmh=args.tmh,
                all_plots=args.all_plots,
                bin_width=args.bin_width,
                y_range=args.y_range,
                tmh_y_range=args.tmh_y_range,
                svg=args.svg,
                force=args.force
            )

        elif args.subparser_name == 'cugo':
            cugo(
                db_path=args.db_path,
                cugo_range=args.cugo_range,
                fasta=args.fasta,
                id_list=args.id_list,
                output_dir=args.output,
                annotation=args.annotation,
                flank_lower=args.flank_lower,
                flank_upper=args.flank_upper,
                top_n=args.top_n,
                threads=args.threads,
                svg=args.svg,
                force=args.force,
                bin_width=args.bin_width,
                y_range=args.y_range,
                tmh_y_range=args.tmh_y_range,
                export=args.export,
                export_dir=args.export_dir,
                homogeneity_threshold=args.homogeneity_threshold,
                homogeneity_window=args.homogeneity_window,
                sequence_frequency_threshold=args.sequence_frequency_threshold
            )

        elif args.subparser_name == 'cugo_select':
            cugo_select(
                context_path=args.context_path,
                position=args.position,
                db_path=args.db_path,
                output=args.output,
                threads=args.threads,
                filter_seqs=args.filter_seqs,
                force=args.force
            )

        elif args.subparser_name == 'filter':
            filter(
                fasta=args.fasta,
                output=args.output,
                db_path=args.db_path,
                threads=args.threads,
                sql=args.sql,
                force=args.force
            )

        ### PARSER FOR CASM FUNCTIONALITIES AND WORKFLOW ###
        elif args.subparser_name == 'casm_plot':
            casm_plot(
                clust_path=args.cluster_path,
                output=args.output,
                metadata=args.metadata,
                db_path=args.db_path,
                force=args.force,
                svg=args.svg,
                show_cluster_numbers=args.show,
            )

        elif args.subparser_name == 'casm':
            casm(
                fasta=args.fasta,
                output=args.output,
                db_path=args.db_path,
                subset=args.subset,
                subset_size=args.subset_size,
                threads=args.threads,
                perplexity=args.perplexity,
                iterations=args.iterations,
                exaggeration=args.exaggeration,
                metadata=args.metadata,
                n_svd_components=args.n_svd_components,
                keep=args.keep,
                force=args.force,
                svg=args.svg,
                show_cluster_numbers=args.show,
                large=args.large
            )

        elif args.subparser_name == 'casm_select':
            casm_select(
                final_embedding_file=args.cluster_path,
                fasta=args.fasta,
                no_cluster=args.no_cluster,
                output=args.output,
                threads=args.threads,
                filter_seqs=args.filter_seqs,
                force=args.force
            )

        elif args.subparser_name == 'list_metadata':
            list_metadata()

        elif args.subparser_name == 'export_fasta':
            export_fasta(
                db_path=args.db_path,
                output=args.output,
                threads=args.threads,
                force=args.force
            )

        elif args.subparser_name == 'export':
            export_by_metadata(
                db_path=args.db_path,
                query=parse_query(args.metadata_query),
                fasta=args.fasta,
                id_list=args.id_list,
                genome_id_list=args.genome_id_list,
                taxonomy_list=args.taxonomy_list,
                taxonomy_column=args.taxonomy_column,
                output=args.output,
                threads=args.threads,
                force=args.force
            )

    except Exception as e:
        logger.error(f"Error executing command: {e}")
        return 1
    return 0

if __name__ == "__main__":
    exit(main())
