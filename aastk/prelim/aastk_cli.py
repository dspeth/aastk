# aastk/cli.py
import argparse
from .pasr import extract_matching_sequences
from .cugo import parse_gff_to_tab, extract_genomic_context

def main():
    parser = argparse.ArgumentParser(description="Amino Acid Sequence Toolkit")

    subparsers = parser.add_subparsers(dest="command", required=True)

    # Subcommand for 'pasr'
    pasr_parser = subparsers.add_parser("pasr", help="PASR: Protein Analysis and Sequence Retrieval")
    pasr_parser.add_argument("blast_tab", help="BLAST/DIAMOND tabular output file")
    pasr_parser.add_argument("read_file", help="FASTA or FASTQ file of sequencing reads")
    pasr_parser.add_argument("out_fasta", help="Output file for matched sequences")

    # Subcommand for 'cugo'
    cugo_parser = subparsers.add_parser("cugo", help="CUGO: Parse GFF files into tab-delimited format")
    cugo_parser.add_argument("-g", "--gff_file", required=True, help="anvi'o generated GFF file to parse")
    cugo_parser.add_argument("-o", "--out_file", required=True, help="Output file for the parsed data")
    cugo_parser.add_argument("-V", "--version", action="store_true", help="Show script version and exit")

    # Subcommand for 'context'
    context_parser = subparsers.add_parser("context", help="Context: Extract genomic context across multiple CUGO files")
    context_parser.add_argument("-i", "--prot_ids", required=True, help="File with list of protein accession numbers")
    context_parser.add_argument("-d", "--cugo_dir", required=True, help="Directory with CUGO files")
    context_parser.add_argument("-t", "--tmhmm_dir", help="OPTIONAL: Directory with TMHMM files")
    context_parser.add_argument("-r", "--cugo_range", type=int, default=0, help="Range of genomic context (default: 0)")
    context_parser.add_argument("-V", "--version", action="store_true", help="Show script version and exit")
    context_parser.add_argument("-o", "--out_file", required=True, help="Output file for genomic context data")

    args = parser.parse_args()

    if args.command == "pasr":
        extract_matching_sequences(args.blast_tab, args.read_file, args.out_fasta)
    elif args.command == "cugo":
    
    elif args.command == "cugo":
        if args.version:
            print("CUGO parser version 0.1")
            print("CUGO: colocated unidirectional gene organization")
        else:
            parse_gff_to_tab(args.gff_file, args.out_file)
    elif args.command == "context":
        if args.version:
            print("Context extractor version 0.1")
        else:
            extract_genomic_context(
                prot_ids=args.prot_ids,
                cugo_dir=args.cugo_dir,
                tmhmm_dir=args.tmhmm_dir,
                cugo_range=args.cugo_range,
                out_file=args.out_file
            )

if __name__ == "__main__":
    main()
