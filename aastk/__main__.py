#!/usr/bin/env python3

from aastk.cli import get_main_parser
from aastk.pasr import *
from aastk.cugo import *

if __name__ == "__main__":
    parser = get_main_parser()
    args = parser.parse_args()

    # If no subcommand is specified, print help
    if not args.subparser_name:
        parser.print_help()

    ### PARSER FOR PASR FUNCTIONALITIES AND WORKFLOW ###
    elif args.subparser_name == 'build':
        db_path = build_protein_db(
            db_dir=args.db,
            protein_name=args.protein_name,
            seed_fasta=args.seed,
            threads=args.threads
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
            chunk=args.chunk
        )

    elif args.subparser_name == 'extract':
        extract_matching_sequences(
            protein_name=args.protein_name,
            blast_tab=args.tabular,
            query_path=args.query,
            output_dir=args.output,
            key_column=args.key_column
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
            score_column=args.score_column
        )

    elif args.subparser_name == 'pasr_plot':
        plot_bsr(
            protein_name=args.protein_name,
            bsr_file=args.bsr,
            output_dir=args.output
        )

    elif args.subparser_name == 'metadata':
        metadata(
            selfmin=args.selfmin,
            selfmax=args.selfmax,
            dbmin=args.dbmin,
            bsr=args.bsr_cutoff,
            output_dir=args.output,
            dataset=args.dataset,
            protein_name=args.protein_name
        )

    elif args.subparser_name == 'subset':
        subset(
            yaml_path=args.yaml,
            matched_fasta=args.matched,
            bsr_table=args.bsr,
            output_dir=args.ouput
        )

    elif args.subparser_name == 'pasr':
        pasr(
            db_dir=args.db,
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
            yaml_path=args.yaml
        )

    ### PARSER FOR CUGO FUNCTIONALITIES AND WORKFLOW ###
    elif args.subparser_name == 'parse':
        parse(
            gff_file_path=args.gff_path,
            output_dir=args.output
        )

    elif args.subparser_name == 'context':
        context(
            protein_ids=args.protein_ids,
            cugo_dir=args.cugo_dir,
            tmhmm_dir=args.tmhmm_dir,
            cugo_range=args.cugo_range,
            output_dir=args.output,
            dataset=args.dataset
        )

    elif args.subparser_name == 'cugo_plot':
        cugo_plot(
            cugo_path=args.cugo_path,
            flank_lower=args.flank_lower,
            flank_upper=args.flank_upper,
            top_n=args.top_n
        )