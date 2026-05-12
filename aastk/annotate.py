# to do: change final output file, so it's based on actual query/pyrodigal FASTA (instead of BSR file)
# to do: plotting like in pasr

from .util import *
from .pasr import *

import logging
import gzip
import pandas as pd
import yaml as yaml_lib

import pyrodigal

logger = logging.getLogger(__name__)

##### functions to be added to util.py #####
def stream_fasta(fasta: str):
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

def clean_fasta_header(header: str):
    """
    keeps the first part of the FASTA header up until the first whitespace
    """
    return header.split()[0]

def file_has_allowed_extension(
    path: str,
    allowed_extensions: tuple[str, ...]):
    """
    returns a boolean whether the file type fits any of the given extensions or not
    """
    return str(path).endswith(allowed_extensions)


##### end of util.py functions #####

def check_fasta_protein(path: str):
    if not file_has_allowed_extension(path, (".faa", ".faa.gz")):
        logger.warning("Protein query input does not have the expected file extension (.faa or .faa.gz)")

    for header, sequence in stream_fasta(path):
        # only checks the first record
        if not header:
            raise ValueError("Protein FASTA contains an empty header")

        allowed_characters = set("ABCDEFGHIKLMNPQRSTVWYZJUX*-")

        invalid_characters = set(sequence.upper()) - allowed_characters

        if invalid_characters:
            raise ValueError("Protein FASTA contains invalid characters")

        return

    raise ValueError("Protein FASTA contains no sequences")

def check_fasta_genome(path: str):
    if not file_has_allowed_extension(path, (".fa", ".fa.gz", ".fna", ".fna.gz")):
        logger.warning("Protein query input does not have the expected file extension (.fa, .fa.gz, .fna or .fna.gz)")

    for header, sequence in stream_fasta(path):
        # only checks the first record
        if not header:
            raise ValueError("Genome FASTA contains an empty header")

        allowed_characters = set("ACGTUIRYKMSWBDHVN-")

        invalid_characters = set(sequence.upper()) - allowed_characters

        if invalid_characters:
            raise ValueError("Genome FASTA contains invalid characters")

        return

    raise ValueError("Genome FASTA contains no sequences")

def copy_protein_query(query: str, output: str, force: bool = False):
    """
    creates a clean, unzipped protein fasta file with a cleaned up header
    """

    check_fasta_protein(query)

    query_name = determine_dataset_name(query, ".", 0)
    output_path = ensure_path(output, f"{query_name}_query.faa", force=force)

    logger.info(f"Input {query} contains amino acid sequences, copying to {output_path}")

    sequence_count = 0

    with open(output_path, "w") as out:
        for header, sequence in stream_fasta(query):
            clean_header = clean_fasta_header(header)
            out.write(f">{clean_header}\n{sequence}\n")
            sequence_count += 1

    logger.info(f"Wrote {sequence_count} protein sequences to {output_path}")

    return output_path

def call_genes_with_pyrodigal(genome: str, output: str, force: bool = False):
    """
    uses pyrodigal to do gene calling and predict amino acid sequences for the input genome/contig
    """
    check_fasta_genome(genome)

    genome_name = determine_dataset_name(genome, ".", 0)
    output_path = ensure_path(output, f"{genome_name}_genes.faa", force=force)
    gff_path = ensure_path(output, f"{genome_name}_genes_info.gff", force=force)

    logger.info(f"Calling genes with Pyrodigal for the input: {genome}")

    gene_finder = pyrodigal.GeneFinder(meta=True)

    contig_count = 0
    protein_count = 0

    with open(output_path, "w") as out_faa, open(gff_path, "w") as out_gff:
        for contig_header, sequence in stream_fasta(genome):
            contig_count += 1
            contig_id = clean_fasta_header(contig_header)

            # pyrodigal analyzes the current sequence and predicts genes
            genes = gene_finder.find_genes(sequence)

            for gene_id, gene in enumerate(genes, start=1):
                protein = gene.translate(include_stop=False)

                if not protein:
                    continue

                protein_id = f"{contig_id}___{gene_id}"

                out_faa.write(f">{protein_id}\n")
                out_faa.write(f"{protein}\n")

                if gene.strand == 1:
                    strand = "+"
                else:
                    strand = "-"

                out_gff.write(
                    f"{contig_id}\t"
                    f"pyrodigal\t"
                    f"CDS\t"
                    f"{gene.begin}\t"
                    f"{gene.end}\t"
                    f".\t"
                    f"{strand}\t"
                    f".\t"
                    f".\t"
                )

                protein_count += 1

    logger.info(f"Predicted {protein_count} protein sequence(s) (out of {contig_count} nucleotide sequence(s))")
    logger.info(f"Saved gene calling results to {output_path}")
    logger.info(f"Saved additional gene calling information to {gff_path}")

    return output_path

def load_yaml_cutoffs(yaml_path: str):
    with open(yaml_path) as f:
        params = yaml_lib.safe_load(f)

    return {
        "bsr_cutoff": float(params["bsr"]),
        "max_score_min": int(params["selfmin"]),
        "max_score_max": int(params["selfmax"]),
        "protein_name": params.get("protein_name"),
    }

def classify_by_yaml_cutoffs(
        max_score: int,
        hit_score: int,
        bsr_cutoff: float,
        max_score_min: int,
        max_score_max: int):

    """
    classify hits based off of the cut-offs in the yaml file
    possible classifications: too_short, correct_length, too_long, below_cutoff
    """

    halfway = max_score_min + ((max_score_max - max_score_min) / 2)

    bsr_line_value = bsr_cutoff * max_score
    halfway_line_value = bsr_cutoff * halfway

    if max_score < max_score_min:
        if hit_score >= bsr_line_value:
            return "too_short"
        return "below_cutoff"

    if max_score_min <= max_score <= halfway:
        if hit_score >= bsr_line_value:
            return "correct_length"
        return "below_cutoff"

    if halfway < max_score <= max_score_max:
        if hit_score >= halfway_line_value:
            return "correct_length"
        return "below_cutoff"

    if max_score > max_score_max:
        if hit_score >= halfway_line_value:
            return "too_long"
        return "below_cutoff"

    raise RuntimeError("Classification of hits failed: max_score is out of bounds")

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
            max_score_min=cutoffs["max_score_min"],
            max_score_max=cutoffs["max_score_max"]
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