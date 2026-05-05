# to do: uniform header for the copy_protein_query and call_genes_with_pyrodigal output FASTAs
# to do: use gff file to save pyrodigal information to keep a clean header in FASTA file
# to do: check filetype (.faa(.gz) or .fa(.gz)) and throw error message if it doesn't fit
# to do: complete workflow only works for protein input so far (maybe fixing pyrodigal output FASTA header will fix that anyways)

from .util import *
from .pasr import *

import logging
import gzip
import pandas as pd
import yaml as yaml_lib

import pyrodigal

logger = logging.getLogger(__name__)

##### functions to be added to util.py #####
def read_fasta_content(fasta: str):
    """
    reads a fasta file entry by entry and returns one header-sequence pair at a time to save memory
    (generator function: "yield" returns values one by one)
    supports multi line FASTA sequences
    """

    filepath = Path(fasta)

    if filepath.suffix == ".gz":
        file_handle = gzip.open(filepath, "rt")
    else:
        file_handle = open(filepath, "r")

    header = None
    sequence_parts = []

    with file_handle as f:
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

                protein_id = f"{contig_id}___{gene_id}"

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

def load_yaml_cutoffs(yaml_path: str):
    with open(yaml_path) as f:
        params = yaml_lib.safe_load(f)

    return {
        "bsr_cutoff": float(params["bsr"]),
        "selfmin": float(params["selfmin"]),
        "selfmax": float(params["selfmax"]),
        "protein_name": params.get("protein_name"),
    }

def classify_by_yaml_cutoffs(
        max_score: float,
        hit_score: float,
        bsr_cutoff: float,
        selfmin: float,
        selfmax: float):

    """
    classify hits based off of the cut-offs in the yaml file
    possible classifications: too_short, correct_length, too_long, below_cutoff
    """

    halfway = selfmin + ((selfmax - selfmin) / 2)

    bsr_line_value = bsr_cutoff * max_score
    halfway_line_value = bsr_cutoff * halfway

    if max_score < selfmin:
        if hit_score >= bsr_line_value:
            return "too_short"
        return "below_cutoff"

    if selfmin <= max_score <= halfway:
        if hit_score >= bsr_line_value:
            return "correct_length"
        return "below_cutoff"

    if halfway < max_score <= selfmax:
        if hit_score >= halfway_line_value:
            return "correct_length"
        return "below_cutoff"

    if max_score > selfmax:
        if hit_score >= halfway_line_value:
            return "too_long"
        return "below_cutoff"

    return "below_cutoff"

def classify_hits(
        bsr_path: str,
        yaml_path: str,
        output: str,
        force: bool = False):

    cutoffs = load_yaml_cutoffs(yaml_path)

    bsr_df = pd.read_csv(bsr_path, sep="\t")

    bsr_df["classification"] = bsr_df.apply(
        lambda row: classify_by_yaml_cutoffs(
            max_score=row["max_score"],
            hit_score=row["score"],
            bsr_cutoff=cutoffs["bsr_cutoff"],
            selfmin=cutoffs["selfmin"],
            selfmax=cutoffs["selfmax"]
        ),
        axis=1
    )

    protein_name = cutoffs["protein_name"]
    output_path = ensure_path(output, f"{protein_name}_annotate.tsv", force=force)

    bsr_df.to_csv(output_path, sep="\t", index=False)

    logger.info(f"Classification complete! Hits saved to: {output_path}")

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

    # get a uniform protein FASTA depending on input
    if query:
        query_path = copy_protein_query(query, output, force=force)

    elif genome:
        query_path = call_genes_with_pyrodigal(genome, output, force=force)

    else:
        raise ValueError("Either query or genome must be provided")

    logger.info(f"Uniform query protein FASTA saved in: {query_path}")

    # build the DIAMOND DB from input db_path
    diamond_db = build(
        seed_fasta=db_path,
        threads=1,
        output=output,
        force=force
    )

    # search query FASTA against the DIAMOND DB
    hits_path, column_info_path = search(
        db_path=diamond_db,
        query_path=query_path,
        threads=1,
        output_dir=output,
        sensitivity="fast",
        force=force
    )

    # calculate max scores for all query proteins
    max_scores_path = max_score(
        extracted=query_path,
        output_dir=output,
        matrix="BLOSUM45",
        force=force
    )

    # calculate bsr for all hits
    bsr_path = bsr(
        blast_tab=hits_path,
        max_scores_path=max_scores_path,
        output_dir=output,
        key_column=0,
        column_info_path=column_info_path,
        force=force
    )

    # classify hits using yaml cutoffs
    annotate_tsv = classify_hits(
        bsr_path=bsr_path,
        yaml_path=yaml,
        output=output,
        force=force
    )

    logger.info("Annotate workflow completed")

    return annotate_tsv