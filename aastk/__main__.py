#!/usr/bin/env python3

import logging
from .log import logger_setup
from aastk.cli import get_main_parser
from aastk.pasr import *
from aastk.cugo import *
from aastk.casm import *
from aastk.database import *
from aastk import __version__, __copyright__, __author__
import sys

def print_help():
    print('''\

              ...::: AASTK v%s :::...

  Workflows:
    pasr (Protein Alignment Score Ratio) -> Generate comprehensive datasets of homologous protein complexes.
                                            (build -> search -> extract -> calculate -> bsr -> pasr_plot)
                                            
    casm  (Clustering Alignment Score Matrix) -> Identify functional heterogeneity within homologous datasets by clustering alignment score matrices.
                                                 (matrix -> cluster -> casm_plot)
                                                 
    cugo (Co-localized Unidirectional Gene Organization) -> Analyze consensus genomic context of selected protein complex clusters.
                                                            (context -> cugo_plot)
 

  Workflow-adjacent methods:
    pasr:
        build -> Build DIAMOND database from seed sequence(s)
        search -> Search DIAMOND reference database for homologous sequences
        extract -> Extract reads that have DIAMOND hits against custom database
        calculate -> Calculate max scores for extracted sequences using BLOSUM matrix
        bsr -> Compute BSR (Blast Score Ratio) using a BLAST tab file and max scores from a TSV.
        pasr_plot -> Plot the Blast Score Ratio of query sequences against the DIAMOND database
        
        Standalone:
        select -> Select target sequences in accordance with metadata cutoffs
        
    casm:
        matrix -> Create alignment matrix for tSNE embedding and DBSCAN clustering
        cluster -> Run tSNE embedding and DBSCAN clustering on input matrix and matrix metadata
        casm_plot -> Plot CASM .tsv output files
        
        Standalone:
        pick -> Pick CASM clusters to generate .faa file for further analysis
        
    cugo:
        context -> Parse context information from CUGO input file
        cugo_plot -> Plot CUGO context
        
        Standalone:
        retrieve -> Retrieve protein IDs for select CUGO position

  Standalone tools:
    meta -> Retrieve metadata from AASTK SQLite database

 


  Use: aastk <command> -h for command specific help; aastk --silent <command> to suppress all console output except errors
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
                svg=args.svg,
                force=args.force,
                update=args.update
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

        elif args.subparser_name == 'retrieve':
            retrieve(
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
                force=args.force,
                svg=args.svg,
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
