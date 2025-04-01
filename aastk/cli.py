#!/usr/bin/env python3

import argparse
from aastk.pasr import build_protein_db
from aastk.pasr import search_protein_db
from aastk.pasr import extract_matching_sequences


def parse_cli():
    parser = argparse.ArgumentParser(description="Amino Acid Sequence Toolkit")

    subparsers = parser.add_subparsers(dest="command", required=True)

    # Subcommand for 'pasr'
    pasr_parser = subparsers.add_parser("pasr", help="PASR: protein alignment score ratio")
    pasr_parser.add_argument("-b", "--blast_tab", required=True, help="BLAST/DIAMOND tabular output file")
    pasr_parser.add_argument("-s", "--seq_file", required=True, help="FASTA or FASTQ file")
    pasr_parser.add_argument("-o", "--out_fasta", required=True, help="Output file for matched sequences")

    db_parser = subparsers.add_parser("build-db", help="Build a DIAMOND database from seed sequence(s)")
    db_parser.add_argument("-p", "--protein_name", required=True, help="Name of the protein for the diamond database")
    db_parser.add_argument("-f", "--seed_fasta", required=True, help="Path to FASTA file containing seed sequence(s)")
    db_parser.add_argument("-t", "--threads", type=int, help="Number of threads to be used by diamond makedb")

    search_parser = subparsers.add_parser("search-db", help="Search DIAMOND reference database for homologs")
    search_parser.add_argument("-d", "--db_path", required=True, help="Path to DIAMOND reference database")
    search_parser.add_argument("-q", "--query_path", required=True, help="Path to query database")
    search_parser.add_argument("-p", "--protein_name", required=True, help="Name of protein of interest")
    search_parser.add_argument("-t", "--threads", type=int, help="Number of threads used by homolog search")

    args = parser.parse_args()

    if args.command == "build-db":
        build_protein_db(args.protein_name, args.seed_fasta, args.threads)

    elif args.command == "search-db":
        search_protein_db(args.db_path, args.query_path, args.threads, args.protein_name)

    elif args.command == "pasr":
        extract_matching_sequences(args.blast_tab, args.seq_file, args.out_fasta)

if __name__ == "__main__":
    parse_cli()