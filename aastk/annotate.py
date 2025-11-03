import logging
import pyrodigal
from Bio import SeqIO
from aastk.util import *



#same logger setup like casm.py
logger = logging.getLogger(__name__)

def gene_calling(
        input_fasta: str,
        output_dir: str,
        detailed_headers: bool = False,
        force: bool = False,
) -> str:
    """
        Predict protein sequences from a nucleotide FASTA using Pyrodigal.

        Args:
            input_fasta (str): Path to the input genome FASTA file.
            output_dir (str): Directory where the predicted protein FASTA will be written.
            detailed_headers (bool): Include coordinates and strand info in FASTA headers.
            force (bool): Overwrite existing output if True.

        Returns:
            str: Path to the generated amino acid FASTA (.faa) file.
        """
    # build dataset name and output path copied from casm
    dataset = determine_dataset_name(input_fasta, splitter=".", part=0)
    output_file = f"{dataset}.proteins.faa"
    output_path = ensure_path(output_dir, output_file, force=force)

    logger.info(f"Starting gene calling with Pyrodigal | input={input_fasta}")

    # initialize and train Pyrodigal model
    finder = pyrodigal.GeneFinder()
    #logger.info("Training Pyrodigal model on complete genome.")
    #finder.train(input_fasta)

    # read genome records (single or multiple contigs)
    records = list(SeqIO.parse(input_fasta, "fasta"))
    total_records = len(records)
    logger.info(f"Found {total_records} sequence record(s) in FASTA.")

    total_genes = 0
    with open(output_path, "w") as out_faa:
        for record in records:
            # predict genes on this sequence
            results = finder.find_genes(str(record.seq))

            for i, gene in enumerate(results):
                if detailed_headers:
                    header = (
                        f">{record.id}_gene{i + 1} "
                        f"{gene.begin}:{gene.end} ({'+' if gene.strand == 1 else '-'})"
                    )
                else:
                    header = f">{record.id}_gene{i + 1}"

                aa_seq = gene.translate()
                out_faa.write(f"{header}\n{aa_seq}\n")

            total_genes += len(results)

    logger.info(
        f"Predicted {total_genes} proteins from {total_records} input sequence(s). "
        f"Output written to {output_path}"
    )

    return output_path


def annotate(
        input_fasta: str,
        output_dir: str,
        detailed_headers: bool = False,
        force: bool = False,
) -> str:
        gene_calling.file = gene_calling(input_fasta, output_dir, detailed_headers, force)
