import csv
import pyrodigal
from aastk.util import *


#same logger setup like casm.py
logger = logging.getLogger(__name__)

def gene_calling(
        input_fasta: str,
        output_dir: str,
        force: bool = False,
) -> str:
    """
        Predict protein sequences from a nucleotide FASTA using Pyrodigal.

        Args:
            input_fasta (str): Path to the input genome FASTA file.
            output_dir (str): Directory where the predicted protein FASTA will be written.
            force (bool): Overwrite existing output if True.

        Returns:
            str: Path to the generated amino acid FASTA (.faa) file.
        """
    # 1) Validate & prepare paths via util
    if determine_file_type(input_fasta) != "fasta":
        raise ValueError(f"Input must be FASTA: {input_fasta}")

    dataset = Path(input_fasta).stem
    output_file = f"{dataset}.proteins.faa"
    output_path = ensure_path(output_dir, output_file, force=force)

    genes_tsv_name = f"{dataset}.genes.tsv"
    genes_tsv_path = ensure_path(output_dir, genes_tsv_name, force=force)

    logger.info(f"Starting gene calling with Pyrodigal | input={input_fasta}")

    # 2) Load sequences with util (no Biopython)
    seqs = read_fasta_to_dict(input_fasta)
    if not seqs:
        raise ValueError(f"No sequences found in {input_fasta}")


    # 3) Initialize and train Pyrodigal model on first sequence (complete-genome assumption)
    first_seq = next(iter(seqs.values()))
    finder = pyrodigal.GeneFinder()
    logger.info("Training Pyrodigal model on complete genome (first record).")
    finder.train(first_seq)

    # 4) Predict and write proteins plus total records
    records = list(seqs.items())
    total_records = len(records)
    logger.info(f"Found {total_records} sequence record(s) in FASTA.")

    total_genes = 0
    with open(output_path, "w") as out_faa, open(genes_tsv_path, "w") as genes_tsv:
        tsv_writer = csv.writer(genes_tsv, delimiter="\t")
        tsv_writer.writerow(["contig_id", "gene_index", "start", "end", "strand", "score"])


        for header, nt_seq in records:
            results = finder.find_genes(nt_seq)
            rec_id = header.split("|")[0]

            for i, gene in enumerate(results, start=1):
                gene_id = f"{rec_id}_{i}"
                strand = "+" if gene.strand == "1" else "-"

                aa_seq = gene.translate()
                aa_length = len(aa_seq)
                out_faa.write(f"{gene_id} | length ={aa_length}AA\n{aa_seq}\n")
                tsv_writer.writerow([rec_id, i, gene.begin, gene.end, strand, gene.score])

            total_genes += 1

    logger.info(
        f"Predicted {total_genes} proteins from {total_records} input sequence(s). "
    )
    logger.info(
        f"Output written to {output_path}"
    )

    return output_path


def annotate(
        input_fasta: str,
        output_dir: str,
        force: bool = False,
) -> str:
    return gene_calling(input_fasta, output_dir, force)
