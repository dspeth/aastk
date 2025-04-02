#!/usr/bin/env python3

from aastk.cli_new import get_main_parser
from aastk.pasr import *

if __name__ == "__main__":
    parser = get_main_parser()
    args = parser.parse_args()

    # If no subcommand is specified, print help
    if not args.subparser_name:
        parser.print_help()
    elif args.subparser_name == 'pasr':
        max_scores = pasr(
            protein_name=args.protein_name,
            seed_fasta=args.seeds,
            query_fasta=args.query,
            matrix_name=args.matrix,
            threads=args.threads,
            target_dir=args.output
        )
        print(f"PASR Results: {max_scores}")