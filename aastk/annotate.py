"""
annotate.py – Step 1: Gene calling with Pyrodigal
Converts a genome FASTA into predicted amino acid sequences (.faa)
"""

import argparse
import pyrodigal


def predict_proteins(input_fasta, output_faa, meta=False, detailed_headers=False):

    # Initialize Pyrodigal gene finder
    finder = pyrodigal.GeneFinder(meta=meta)

    # Train if not in meta mode (required for complete genomes)
    if not meta:
        finder.train(input_fasta)

    # Run gene prediction
    result = finder.find_genes(input_fasta)

    # Write predicted proteins (AA sequences) to FASTA
    with open(output_faa, "w") as out:
        for i, gene in enumerate(result):
            # Header formatting
            if detailed_headers:
                header = (
                    f">{gene.id}_gene{i+1} "
                    f"{gene.begin}:{gene.end} ({'+' if gene.strand == 1 else '-'})"
                )
            else:
                header = f">{gene.id}_gene{i+1}"

            seq = gene.translation
            out.write(f"{header}\n{seq}\n")

    print(
        f"Predicted {len(result)} proteins written to '{output_faa}' "
        f"(meta={meta}, detailed_headers={detailed_headers})"
    )


def main():
    """Handle command-line interface (CLI) arguments."""
    parser = argparse.ArgumentParser(
        description="Predict protein-coding genes from a genome FASTA using Pyrodigal."
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input genome FASTA file."
    )

    parser.add_argument(
        "-o", "--output",
        default="predicted_proteins.faa",
        help="Output protein FASTA file (default: predicted_proteins.faa)."
    )

    parser.add_argument(
        "--meta",
        action="store_true",
        help="Use metagenomic mode (meta=True). Default is False."
    )

    parser.add_argument(
        "--detailed-headers",
        action="store_true",
        help="Include genomic coordinates and strand info in FASTA headers."
    )

    args = parser.parse_args()

    # Call the main prediction function
    predict_proteins(
        input_fasta=args.input,
        output_faa=args.output,
        meta=args.meta,
        detailed_headers=args.detailed_headers,
    )


if __name__ == "__main__":
    main()
