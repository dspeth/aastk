#!/usr/bin/env python3

import argparse
from aastk.pasr import build_protein_db
from aastk.pasr import search_protein_db
from aastk.pasr import extract_matching_sequences
from aastk.pasr import calculate_max_scores

# look at gtdbtk for cleaner parsing
def parse_cli():
    parser = argparse.ArgumentParser(description="Amino Acid Sequence Toolkit")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # top level command for 'pasr'
    pasr_parser = subparsers.add_parser("pasr", help="PASR: protein alignment score ratio")
    pasr_subparsers = pasr_parser.add_subparsers(dest="subcommand", required=True)

    # subparser for pasr build-db
    build_db_parser = pasr_subparsers.add_parser("build-db", help="Build a DIAMOND database")
    build_db_parser.add_argument("-p", "--protein_name", required=True, help="Protein name for the database")
    build_db_parser.add_argument("-f", "--seed_fasta", required=True, help="Seed FASTA file")
    build_db_parser.add_argument("-t", "--threads", type=int, help="Threads for makedb")

    # subparser for pasr search-db
    search_db_parser = pasr_subparsers.add_parser("search-db", help="Search DIAMOND database")
    search_db_parser.add_argument("-d", "--db_path", required=True, help="Path to DIAMOND reference database")
    search_db_parser.add_argument("-q", "--query_path", required=True, help="Query FASTA file")
    search_db_parser.add_argument("-p", "--protein_name", required=True, help="Protein of interest")
    search_db_parser.add_argument("-t", "--threads", type=int, help="Threads for search")

    # subparser for pasr extract
    extract_parser = pasr_subparsers.add_parser("extract", help="Extract matched sequences")
    extract_parser.add_argument("-b", "--blast_tab", required=True, help="BLAST/DIAMOND tabular output file")
    extract_parser.add_argument("-s", "--seq_file", required=True, help="FASTA or FASTQ file")
    extract_parser.add_argument("-o", "--out_fasta", required=True, help="Output file for matched sequences")

    # subparser for pasr max - score
    max_score_parser = pasr_subparsers.add_parser("max-score", help="Calculate maximum BLOSUM alignment score for sequences")
    max_score_parser.add_argument("-f", "--fasta", required=True, help="Protein FASTA file")
    max_score_parser.add_argument("-m", "--matrix", required=True, choices=["BLOSUM45", "BLOSUM62"], help="BLOSUM matrix to use (BLOSUM45 or BLOSUM62)")
    max_score_parser.add_argument("-o", "--output", required=True, help="Output file for scores")

    # subparser for pasr all
    all_parser = pasr_subparsers.add_parser("all", help="Run all steps: build-db, search-db, extract, max-score")
    all_parser.add_argument("-p", "--protein_name", required=True, help="Protein name for the database")
    all_parser.add_argument("-f", "--seed_fasta", required=True, help="Seed FASTA file")
    all_parser.add_argument("-q", "--query_path", required=True, help="Query FASTA file")
    all_parser.add_argument("-t", "--threads", type=int, required=True, help="Threads for all steps")
    all_parser.add_argument("-m", "--matrix", required=True, choices=["BLOSUM45", "BLOSUM62"], help="BLOSUM matrix to use (BLOSUM45 or BLOSUM62)")
    all_parser.add_argument("-s", "--matched", required=True, help="Path to matched sequences FASTA")
    all_parser.add_argument("-o", "--output", required=True, help="Output file for max score results")

    args = parser.parse_args()

    # check if pasr functionalities should be run
    if args.command == "pasr":
        if args.subcommand == "build-db":
            build_protein_db(args.protein_name, args.seed_fasta, args.threads)

        elif args.subcommand == "search-db":
            search_protein_db(args.db_path, args.query_path, args.threads, args.protein_name)

        elif args.subcommand == "extract":
            extract_matching_sequences(args.blast_tab, args.seq_file, args.out_fasta)

        elif args.subcommand == "max-score":
            calculate_max_scores(args.fasta, args.matrix, args.output)

        elif args.subcommand == "all":
            db_path = build_protein_db(args.protein_name, args.seed_fasta, args.threads)
            blast_tab = search_protein_db(db_path, args.query_path, args.threads, args.protein_name)
            fasta = extract_matching_sequences(blast_tab, args.query_path, args.matched)
            calculate_max_scores(fasta, args.matrix, args.output)

if __name__ == "__main__":
    parse_cli()