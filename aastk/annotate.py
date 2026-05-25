# for directory input: handle multiple hits in classification_output?

from .util import *
from .pasr import *

import logging
import gzip
import pandas as pd
import yaml as yaml_lib

import pyrodigal

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

logger = logging.getLogger(__name__)

##### functions to maybe be added to util.py #####
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
    for header, sequence in stream_fasta(path):
        # only checks the first record
        if not header:
            raise ValueError("Protein FASTA contains an empty header")

        allowed_characters = set("ABCDEFGHIKLMNPQRSTVWYZJUX*")

        invalid_characters = set(sequence.upper()) - allowed_characters

        if invalid_characters:
            raise ValueError("Protein FASTA contains invalid characters")

        return

    raise ValueError("Protein FASTA contains no sequences")

def check_fasta_genome(path: str, acgt_threshold: float = 0.8):
    for header, sequence in stream_fasta(path):
        # only checks the first record
        if not header:
            raise ValueError("Genome FASTA contains an empty header")

        allowed_characters = set("ACGTUIRYKMSWBDHVN")

        invalid_characters = set(sequence.upper()) - allowed_characters

        if invalid_characters:
            raise ValueError("Genome FASTA contains invalid characters")

        standard_bases = set("ACGT")
        standard_base_count = sum(1 for char in sequence if char in standard_bases)
        standard_base_ratio = standard_base_count / len(sequence)

        if standard_base_ratio < acgt_threshold:
            raise ValueError(f"Genome FASTA contains too few standard nucleic bases, expected ratio >{acgt_threshold}")

        return

    raise ValueError("Genome FASTA contains no sequences")

def clean_protein_query(query: str, output: str, force: bool = False):
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

            genes.write_gff(
                out_gff,
                sequence_id=contig_id,
                header=(contig_count == 1),
                include_translation_table=True
            )

            for gene_id, gene in enumerate(genes, start=1):
                protein = gene.translate(include_stop=False)

                if not protein:
                    continue

                protein_id = f"{contig_id}___{gene_id}"

                out_faa.write(f">{protein_id}\n")
                out_faa.write(f"{protein}\n")

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
        "max_score_lower": int(params["selfmin"]),
        "max_score_upper": int(params["selfmax"]),
        "protein_name": params.get("protein_name"),
    }

def classify_by_yaml_cutoffs(
        max_score: int,
        hit_score: int,
        bsr_cutoff: float,
        max_score_lower: int,
        max_score_upper: int):

    """
    classify hits based off of the cut-offs in the yaml file
    possible classifications: too_short, correct_length, too_long, below_cutoff
    """

    halfway = max_score_lower + ((max_score_upper - max_score_lower) / 2)

    bsr_line_value = bsr_cutoff * max_score
    halfway_line_value = bsr_cutoff * halfway

    if max_score < max_score_lower:
        if hit_score >= bsr_line_value:
            return "too_short"
        return "below_cutoff"

    if max_score_lower <= max_score <= halfway:
        if hit_score >= bsr_line_value:
            return "correct_length"
        return "below_cutoff"

    if halfway < max_score <= max_score_upper:
        if hit_score >= halfway_line_value:
            return "correct_length"
        return "below_cutoff"

    if max_score > max_score_upper:
        if hit_score >= halfway_line_value:
            return "too_long"
        return "below_cutoff"

    raise RuntimeError("Classification of hits failed: max_score is out of bounds")

def classification_output(
        query_path: str,
        bsr_path: str,
        yaml_path: str,
        output: str,
        force: bool = False):

    """
    classifies the hits based off of classify_by_yaml_cutoffs,
    creates an output file that includes hits and non-hits
    """

    cutoffs = load_yaml_cutoffs(yaml_path)
    protein_name = cutoffs["protein_name"]

    query_df = create_output_rows(query_path)
    bsr_df = pd.read_csv(bsr_path, sep="\t")

    query_id_column = bsr_df.columns[0]

    bsr_df["classification"] = bsr_df.apply(
        lambda row: classify_by_yaml_cutoffs(
            max_score=row["max_score"],
            hit_score=row["score"],
            bsr_cutoff=cutoffs["bsr_cutoff"],
            max_score_lower=cutoffs["max_score_lower"],
            max_score_upper=cutoffs["max_score_upper"]
        ),
        axis=1
    )

    bsr_df = bsr_df.rename(columns={query_id_column: "prot_ID"})

    result_df = query_df.merge(
        bsr_df,
        on="prot_ID",
        how="left",
    )

    result_df["classification"] = result_df["classification"].fillna("no_hit")

    result_df["annotation"] = result_df["classification"].apply(
        lambda classification:
            protein_name
            if classification not in ["below_cutoff", "no_hit"]
            else pd.NA
    )

    first_columns = [
        "prot_ID",
        "annotation",
        "classification",
        "sequence_length"
    ]

    additional_columns = []

    for column in result_df.columns:
        if column not in first_columns:
            additional_columns.append(column)

    result_df = result_df[first_columns + additional_columns]

    output_path = ensure_path(output, f"{protein_name}_annotate.tsv", force=force)

    result_df.to_csv(output_path, sep="\t", index=False)

    logger.info(f"Classification complete! Hits saved to: {output_path}")

    return output_path

def create_output_rows(query_path: str):
    """
    creates one output row per protein sequence in the query fasta (hit and non-hit)
    """

    rows = []

    for header, sequence in stream_fasta(query_path):
        prot_id = clean_fasta_header(header)

        rows.append({
            "prot_ID": prot_id,
            "sequence_length": len(sequence)
        })

    return pd.DataFrame(rows)

def annotate_plot(
        bsr_path: str,
        yaml_path: str,
        output: str,
        svg: bool = False,
        force: bool = False):
    """
    creates a scatter plot based on the bsr file,
    output default is png, but can be set to svg
    """

    cutoffs = load_yaml_cutoffs(yaml_path)

    protein_name = cutoffs["protein_name"]
    bsr_cutoff = cutoffs["bsr_cutoff"]
    max_score_lower = cutoffs["max_score_lower"]
    max_score_upper = cutoffs["max_score_upper"]

    halfway = max_score_lower + ((max_score_upper - max_score_lower) / 2)
    halfway_line_value = bsr_cutoff * halfway

    bsr_df = pd.read_csv(bsr_path, sep="\t")

    if bsr_df.empty:
        logger.warning("No hits found in BSR table. Skipping annotate plot.")
        return None

    required_columns = [
        "max_score",
        "score",
        "pident",
        "BSR",
    ]

    for column in required_columns:
        if column not in bsr_df.columns:
            raise ValueError(f"Required column '{column}' not found in BSR table")

    if svg:
        output_path = ensure_path(output, f"{protein_name}_annotate_plot.svg", force=force)
    else:
        output_path = ensure_path(output, f"{protein_name}_annotate_plot.png", force=force)

    logger.info(f"Creating annotate plot for {protein_name}")

    bin_width = 10

    score_max = bsr_df["score"].max()
    max_score_max = bsr_df["max_score"].max()

    x_max = max(max_score_max * 1.2, max_score_upper * 1.2)
    y_max = max(score_max * 1.2, halfway_line_value * 1.5)

    xlim = (0, x_max)
    ylim = (0, y_max)

    x_bins = np.arange(xlim[0], xlim[1] + bin_width, bin_width)
    y_bins = np.arange(ylim[0], ylim[1] + bin_width, bin_width)

    fig, axs = plt.subplot_mosaic(
        [
            ["histx", "."],
            ["scatter", "histy"],
        ],
        figsize=(8, 8),
        width_ratios=(4, 1),
        height_ratios=(1, 4),
        layout="constrained",
    )

    scatter = axs["scatter"].scatter(
        bsr_df["max_score"],
        bsr_df["score"],
        c=bsr_df["pident"],
        cmap="viridis",
        edgecolor="k",
        alpha=0.5,
        vmin=0,
        vmax=100,
    )

    axs["scatter"].set_xlabel("Calculated maximum score")
    axs["scatter"].set_ylabel("Alignment score to seed set")
    axs["scatter"].set_xlim(xlim)
    axs["scatter"].set_ylim(ylim)

    cb_ax = inset_axes(
        axs["scatter"],
        width="5%",
        height="30%",
        loc="upper left",
        borderpad=1,
    )

    cbar = fig.colorbar(scatter, cax=cb_ax)
    cbar.set_label("% sequence identity")

    x_hist, x_edges = np.histogram(bsr_df["max_score"], bins=x_bins)

    axs["histx"].bar(
        x_edges[:-1],
        x_hist,
        width=bin_width,
        align="edge",
        color="black",
    )

    axs["histx"].set_xlim(xlim)
    axs["histx"].set_ylabel("Counts")
    axs["histx"].tick_params(labelbottom=False)
    axs["histx"].set_title(f"Annotate classification plot for {protein_name}")

    y_hist, y_edges = np.histogram(bsr_df["score"], bins=y_bins)

    axs["histy"].barh(
        y_edges[:-1],
        y_hist,
        height=bin_width,
        align="edge",
        color="black",
    )

    axs["histy"].set_ylim(ylim)
    axs["histy"].set_xlabel("Counts")
    axs["histy"].tick_params(labelleft=False)

    axs["scatter"].axvline(
        max_score_lower,
        color="black",
        linestyle="--",
        linewidth=1.0,
    )

    axs["scatter"].axvline(
        max_score_upper,
        color="black",
        linestyle="--",
        linewidth=1.0,
    )

    x_values_left = np.linspace(0, halfway, 500)
    y_values_left = bsr_cutoff * x_values_left

    axs["scatter"].plot(
        x_values_left,
        y_values_left,
        color="black",
        linestyle="--",
        linewidth=1.0,
    )

    x_values_right = np.linspace(halfway, xlim[1], 500)
    y_values_right = np.full_like(x_values_right, halfway_line_value)

    axs["scatter"].plot(
        x_values_right,
        y_values_right,
        color="black",
        linestyle="--",
        linewidth=1.0,
    )

    fig.savefig(output_path, dpi=300)
    plt.close(fig)

    logger.info(f"Annotate plot saved to: {output_path}")

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
        query_path = clean_protein_query(query, output, force=force)

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
    annotate_tsv = classification_output(
        query_path=query_path,
        bsr_path=bsr_path,
        yaml_path=yaml,
        output=output,
        force=force
    )

    annotate_plot_path = annotate_plot(
        bsr_path=bsr_path,
        yaml_path=yaml,
        output=output,
        force=force
    )

    logger.info("Annotate workflow completed")

    return annotate_tsv