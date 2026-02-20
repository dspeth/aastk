#!/usr/bin/env python3

import logging
from .log import logger_setup
from aastk.cli import get_main_parser
from aastk.pasr import *
from aastk.cugo import *
from aastk.casm import *
from aastk.database import *
from aastk.rasr import *
from aastk import __version__, __copyright__, __author__
import sys
from aastk.rasr_multiple import *

def print_help():
    print('''\

            ...::: AASTK v%s :::...

  Usage:
    aastk <command> <arguments>    to run the tools or commands
    aastk <command> -h             for command specific help
    aastk --silent <command>       to suppress all console output except errors

  Main Tools:
    pasr         Generate or update a comprehensive dataset of homologous proteins
                    Runs the subcommands: build, search, get_hit_seqs, max_score, bsr, and pasr_plot
    casm         Cluster and visualize a complex dataset of proteins (eg. a superfamily)
                    Runs the subcommands: matrix, cluster, and casm_plot
    cugo         Retrieve, calculate, and visualise consensus genomic context of protein data sets
                    Runs the subcommands: context, cugo_plot
    meta         Retrieve protein metadata from AASTK SQLite database

  Helper tools:
    pasr_select  Select target sequences from pasr run based on bsr and score cutoffs
    casm_select  Select cluster(s) from casm analysis, and retrieve sequences
    cugo_select  Select genomic posiiton from cugo analysis and retrieve sequences
    filter       Filter sequence dataset to remove non-homologous seqeunces

  Subcommands:
    build        Build DIAMOND database from seed sequence(s)
    search       Search query sequences against DIAMOND database
    get_hit_seqs Extract sequences that have DIAMOND hits against custom database
    max_score    Calculate max scores for extracted sequences using BLOSUM matrix
    bsr          Compute BSR (Blast Score Ratio) using a BLAST tab file and max scores from a TSV
    pasr_plot    Scatterplot with max score on the x-axis and score against the seed db on y-axis
    matrix       Create alignment score matrix for tSNE embedding and DBSCAN clustering
    cluster      Run tSNE embedding and DBSCAN clustering on input matrix
    casm_plot    Scatterplot of sequences with tSNE coordinates as axes
    context      Parse context information from AASTK SQL database
    cugo_plot    Consensus genomic context plot of annotation, length, and transmembrane segments
    ''' % __version__)


def main():
    parser = get_main_parser()

    if len(sys.argv) == 1:
        print_help()
        sys.exit(0)
    elif sys.argv[1] in {'-v', '--v', '-version', '--version'}:
        print(f"AASTK: version {__version__} {__copyright__} {__author__}")
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

        elif args.subparser_name == 'get_hit_seqs':
            get_hit_seqs(
                blast_tab=args.tabular,
                query_path=args.query,
                output_dir=args.output,
                db_path=args.db_path,
                key_column=args.key_column,
                sql=args.sql,
                force=args.force
            )

        elif args.subparser_name == 'max_score':
            max_score(
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
                svg=args.svg,
                force=args.force,
                update=args.update
            )

        elif args.subparser_name == 'pasr_select':
            pasr_select(
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
                keep=args.keep,
                svg=args.svg,
                force=args.force
            )

        ### PARSER FOR CUGO FUNCTIONALITIES AND WORKFLOW ###
        elif args.subparser_name == 'database':
            database(
                cog_gff_tar_path=args.cog_gff,
                kegg_gff_tar_path=args.kegg_gff,
                pfam_gff_tar_path=args.pfam_gff,
                tmhmm_tar_path=args.tmhmm_dir,
                protein_fasta_path=args.all_proteins,
                taxonomy_path=args.taxonomy_path,
                culture_collection_path=args.culture_collection_path,
                high_level_environment_path=args.high_level_environment_path,
                low_level_environment_path=args.low_level_environment_path,
                output_dir=args.output,
                globdb_version=args.globdb_version
            )

        elif args.subparser_name == 'meta':
            meta(
                db_path=args.db_path,
                fasta=args.fasta,
                output=args.output,
                threads=args.threads,
                include_annotation=args.include_annotation,
                include_taxonomy=args.include_taxonomy,
                all_metadata=args.all_metadata,
                force=args.force
            )

        elif args.subparser_name == 'context':
            context(
                fasta=args.fasta,
                id_list=args.id_list,
                cugo_path=args.cugo_path,
                cugo_range=args.cugo_range,
                output_dir=args.output,
                threads=args.threads,
                force=args.force
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
                all_plots=args.all_plots,
                bin_width=args.bin_width,
                y_range=args.y_range,
                tmh_y_range=args.tmh_y_range,
                svg=args.svg,
                force=args.force
            )

        elif args.subparser_name == 'cugo':
            cugo(
                cugo_path=args.cugo_path,
                cugo_range=args.cugo_range,
                fasta=args.fasta,
                id_list=args.id_list,
                output_dir=args.output,
                flank_lower=args.flank_lower,
                flank_upper=args.flank_upper,
                top_n=args.top_n,
                threads=args.threads,
                svg=args.svg,
                force=args.force,
                bin_width=args.bin_width,
                y_range=args.y_range,
                tmh_y_range=args.tmh_y_range
            )

        elif args.subparser_name == 'cugo_select':
            cugo_select(
                context_path=args.context_path,
                position=args.position,
                db_path=args.db_path,
                output=args.output,
                force=args.force
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
                force=args.force,
            )


        elif args.subparser_name == 'casm_plot':
            casm_plot(
                early_clust_path=args.early_clust,
                full_clust_path=args.full_clust,
                output=args.output,
                db_path=args.db_path,
                metadata_protein=args.metadata_protein,
                metadata_genome=args.metadata_genome,
                force=args.force,
                svg=args.svg,
                show_cluster_numbers=args.show
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
                metadata_protein=args.metadata_protein,
                metadata_genome=args.metadata_genome,
                keep=args.keep,
                force=args.force,
                svg=args.svg,
                show_cluster_numbers=args.show
            )

        elif args.subparser_name == 'casm_select':
            casm_select(
                final_embedding_file=args.full_clust,
                fasta=args.fasta,
                no_cluster=args.no_cluster,
                output=args.output,
                force=args.force
            )

        ### PARSER FOR RASR FUNCTIONALITIES AND WORKFLOW ###
        elif args.subparser_name == 'rasr':
            rasr(
                query=args.query,
                gene_db_fasta=args.seed,
                outgrp_db=args.outgrp_db,
                output_dir=args.output,
                threads=args.threads,
                sensitivity=args.sensitivity,
                block=args.block,
                chunk=args.chunk,
                keep=args.keep,
                force=args.force
            )
        
        elif args.subparser_name == 'rasr_plot':
            rasr_plot(
                bsr_file=args.bsr,
                output_dir=args.output,
                force=args.force,
                bsr_cutoff=args.bsr_cutoff,
                dbmin=args.dbmin
            )

        elif args.subparser_name == 'rasr_select':
            rasr_select(
                score_cutoff=args.dbmin,
                bsr_cutoff=args.bsr_cutoff,
                matched_fastq=args.matched,
                bsr_file=args.bsr,
                output_dir=args.output,
                force=args.force
            )

        ### PARSER FOR RASR_MULTIPLE (Multi-dataset/Multi-gene workflow) ###
        elif args.subparser_name == 'rasr_multiple':
            query_input = args.query_dir if args.query_dir is not None else args.query
            seed_input = args.seed_dir if args.seed_dir is not None else args.seed
            rasr_multiple(
                query=query_input,
                seed_db=seed_input,
                outgrp_db=args.outgrp_db,
                output_dir=args.output,
                sensitivity=args.sensitivity,
                block=args.block,
                chunk=args.chunk,
                threads=args.threads,
                keep=args.keep,
                force=args.force,
                bsr_cutoff=args.bsr_cutoff,
                dbmin=args.dbmin,
                use_existing_merged=args.use_existing_merged
            )


    except Exception as e:
        logger.error(f"Error executing command: {e}")
        return 1
    return 0

if __name__ == "__main__":
    exit(main())
