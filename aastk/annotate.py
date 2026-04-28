# to do: uniform header for the copy_protein_query and call_genes_with_pyrodigal output FASTAs
# to do: error message if -q xyz.fa.gz or --genome xyz.faa.gz?

from .util import *

import logging
import gzip

import pyrodigal

logger = logging.getLogger(__name__)

##### functions to be added to util.py #####

def open_gzip(path: str, mode: str = "rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)

def read_fasta_content(fasta: str):
    """
    reads a fasta file entry by entry and returns one header-sequence pair at a time to save memory
    (generator function: "yield" returns values one by one)
    """

    header = None
    sequence_parts = []

    with open_gzip(fasta, "rt") as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence_parts)

                header = line[1:]
                sequence_parts = []

            else:
                sequence_parts.append(line)
    if header is not None:
        yield header, "".join(sequence_parts)

##### end of util.py functions #####

def copy_protein_query(query: str, output: str, force: bool = False):
    """
    creates a clean, unzipped protein fasta file
    """

    query_name = determine_dataset_name(query, ".", 0)
    output_path = ensure_path(output, f"{query_name}_query.faa", force=force)

    logger.info(f"Input {query} contains amino acid sequences, copying to {output_path}")

    sequence_count = 0

    with open(output_path, "w") as out:
        for header, sequence in read_fasta_content(query):
            out.write(f">{header}\n{sequence}\n")
            sequence_count += 1

    logger.info(f"Wrote {sequence_count} protein sequences to {output_path}")

    return output_path

def call_genes_with_pyrodigal(genome: str, output: str, force: bool = False):
    """
    uses pyrodigal to do gene calling and predict amino acid sequences for the input genome/contig
    """
    genome_name = determine_dataset_name(genome, ".", 0)
    output_path = ensure_path(output, f"{genome_name}_genes.faa", force=force)

    logger.info(f"Calling genes with Pyrodigal for the input: {genome}")

    # creating the pyrodigal gene caller
    # meta means no model is trained on the input first (useful for metagenomes/contigs) -> useable as default? entscheidungsfindung?
    gene_finder = pyrodigal.GeneFinder(meta=True)

    contig_count = 0
    protein_count = 0

    with open(output_path, "w") as out:
        for contig_id, sequence in read_fasta_content(genome):
            contig_count += 1

            # pyrodigal analyzes the current sequence and predicts genes
            genes = gene_finder.find_genes(sequence)

            for gene_id, gene in enumerate(genes, start=1):
                protein = gene.translate(include_stop=False)

                if not protein:
                    continue

                protein_id = f"{contig_id}_{gene_id}"

                out.write(
                    f">{protein_id} "
                    f"begin={gene.begin} "
                    f"end={gene.end} "
                    f"strand={gene.strand}\n"
                )
                out.write(f"{protein}\n")

                protein_count += 1

    logger.info(f"Predicted {protein_count} protein sequence(s) (out of {contig_count} nucleotide sequence(s))")

    return output_path


def annotate(
        db_path: str,
        yaml: str,
        output: str,
        query: str | None = None,
        genome: str | None = None,
        keep: bool = False,
        force: bool = False):

    logger.info("Starting annotate workflow")

    if query:
        query_path = copy_protein_query(query, output, force=force)

    else:
        query_path = call_genes_with_pyrodigal(genome, output, force=force)

    logger.info(f"Uniform query protein FASTA saved in: {query_path}")

    return query_path