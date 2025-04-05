#!/usr/bin/env python3

from aastk.cli_new import get_main_parser
from aastk.pasr import *

if __name__ == "__main__":
    parser = get_main_parser()
    args = parser.parse_args()

    # If no subcommand is specified, print help
    if not args.subparser_name:
        parser.print_help()

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
            sensitivity=args.sensitivity
        )

    elif args.subparser_name == 'extract':
        extract_matching_sequences(
            blast_tab=args.tabular,
            query_path=args.query,
            output_dir=args.output,
            key_column=args.key_column
        )

    elif args.subparser_name == 'calculate':
        calculate_max_scores(
            extracted=args.extracted,
            matrix=args.matrix,
            output_dir=args.output
        )

    elif args.subparser_name == 'pasr':
        max_scores = pasr(
            db_dir=args.db,
            protein_name=args.protein_name,
            seed_fasta=args.seed,
            query_fasta=args.query,
            matrix_name=args.matrix,
            threads=args.threads,
            output_dir=args.output
        )
        print(f"PASR Results: {max_scores}")