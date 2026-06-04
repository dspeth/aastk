# does workflow only provide BEST hit?
# ADJUST ANNOTATE FUNCTION

from .util import *
from .pasr import *

import logging
import gzip
import pandas as pd
import yaml as yaml_lib

import pyrodigal

import matplotlib.pyplot as plt
import numpy as np

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

def annotate_plot_helper(
        bsr_df: pd.DataFrame,
        yaml_path: str,
        output: str,
        marker: str,
        svg: bool = False,
        force: bool = False):
    """
    creates a scatter plot for the given marker,
    output default is png, but can be set to svg
    """

    cutoffs = load_yaml_cutoffs(yaml_path)

    protein_name = cutoffs["protein_name"]
    bsr_cutoff = cutoffs["bsr_cutoff"]
    max_score_lower = cutoffs["max_score_lower"]
    max_score_upper = cutoffs["max_score_upper"]

    halfway = max_score_lower + ((max_score_upper - max_score_lower) / 2)
    halfway_line_value = bsr_cutoff * halfway

    if bsr_df.empty:
        logger.warning(f"No hits found for marker {marker}. Skipping annotate plot.")
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

    score_max = bsr_df["score"].max()
    max_score_max = bsr_df["max_score"].max()

    x_max = max(max_score_max * 1.2, max_score_upper * 3)
    y_max = max(score_max * 1.2, halfway_line_value * 1.5)

    xlim = (0, x_max)
    ylim = (0, y_max)

    fig, ax = plt.subplots(
        figsize=(8, 6),
        layout="constrained"
    )

    scatter = ax.scatter(
        bsr_df["max_score"],
        bsr_df["score"],
        c=bsr_df["pident"],
        cmap="viridis",
        edgecolor="k",
        alpha=0.5,
        vmin=0,
        vmax=100,
    )

    ax.set_title(f"Annotate classification plot for {protein_name} ({marker})")
    ax.set_xlabel("Calculated maximum score")
    ax.set_ylabel("Alignment score to seed set")
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    cbar = fig.colorbar(scatter, ax=ax)
    cbar.set_label("% sequence identity")

    ax.axvline(
        max_score_lower,
        color="black",
        linestyle="--",
        linewidth=1.0,
    )

    ax.axvline(
        max_score_upper,
        color="black",
        linestyle="--",
        linewidth=1.0,
    )

    x_values_left = np.linspace(0, halfway, 500)
    y_values_left = bsr_cutoff * x_values_left

    ax.plot(
        x_values_left,
        y_values_left,
        color="black",
        linestyle="--",
        linewidth=1.0,
    )

    x_values_right = np.linspace(halfway, xlim[1], 500)
    y_values_right = np.full_like(x_values_right, halfway_line_value)

    ax.plot(
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

def marker_name_from_file(path: str):
    path = Path(path)

    if path.suffix == ".gz":
        return Path(path.stem).stem

    return Path(path).stem

def collect_db_fastas(db_path: str):
    """
    collect all fasta files from db_path
    expected folder structure:
    db_path/
        marker1/
            marker1.faa
            marker1.yaml
        marker2/
            marker2.faa
            marker2.yaml
    marker name taken from subdirectory name
    """
    path = Path(db_path)

    if path.is_file():
        check_fasta_protein(str(path))
        marker_name = marker_name_from_file(str(path))
        return {marker_name: str(path)}

    if not path.is_dir():
        raise ValueError(f"Database path does not exist: {db_path}")

    fasta_files = {}

    for marker_dir in sorted(path.iterdir()):
        if not marker_dir.is_dir():
            continue

        marker_name = marker_dir.name
        marker_fastas = []

        for file in sorted(marker_dir.iterdir()):
            if not file.is_file():
                continue

            try:
                check_fasta_protein(str(file))
                marker_fastas.append(file)

            except ValueError:
                continue

        if not marker_fastas:
            logger.warning(f"No valid FASTA file found for marker: {marker_name}")
            continue

        if len(marker_fastas) > 1:
            raise ValueError(f"Multiple FASTA files found for marker: {marker_name}")

        fasta_files[marker_name] = str(marker_fastas[0])

    if not fasta_files:
        raise ValueError(f"No valid FASTA files found in database directory: {db_path}")

    return fasta_files

def collect_yaml_files(db_path: str):
    """
    collects all yaml files from db_path
    expected folder structure:
    db_path/
        marker1/
            marker1.faa
            marker1.yaml
        marker2/
            marker2.faa
            marker2.yaml
    marker name taken from subdirectory name
    """

    path = Path(db_path)

    if path.is_file():
        load_yaml_cutoffs(str(path))
        marker_name = marker_name_from_file(str(path))
        return {marker_name: str(path)}

    if not path.is_dir():
        raise ValueError(f"Database path does not exist: {db_path}")

    yaml_files = {}

    for marker_dir in sorted(path.iterdir()):
        if not marker_dir.is_dir():
            continue

        marker_name = marker_dir.name
        marker_yamls = []

        for file in sorted(marker_dir.iterdir()):
            if not file.is_file():
                continue

            try:
                load_yaml_cutoffs(str(file))
                marker_yamls.append(file)

            except (ValueError, KeyError, TypeError, yaml_lib.YAMLError):
                continue

        if not marker_yamls:
            logger.warning(f"No valid YAML files found for marker: {marker_name}")
            continue

        if len(marker_yamls) > 1:
            raise ValueError(f"Multiple YAML files found for marker: {marker_name}")

        yaml_files[marker_name] = str(marker_yamls[0])

    if not yaml_files:
        raise ValueError(f"No valid YAML files found in database directory: {db_path}")

    return yaml_files

def check_db_yaml_matching(db_markers: set[str], yaml_markers: set[str]):
    """
    checks if every db marker has a matching yaml file and vice versa
    """

    missing_yamls = db_markers - yaml_markers
    missing_dbs = yaml_markers - db_markers

    if missing_yamls:
        missing_yamls_print = ", ".join(sorted(missing_yamls))

        raise ValueError(
            f"Missing YAML file(s) for: {missing_yamls_print}"
        )

    if missing_dbs:
        missing_dbs_print = ", ".join(sorted(missing_dbs))

        raise ValueError(
            f"YAML file(s) without matching db FASTA: {missing_dbs_print}"
        )

def concatenate_db_input(db_fastas: dict[str, str], output: str, force: bool = False):
    """
    concatenates all db input files into a single fasta (necessary for DIAMOND search)
    keeps the headers/IDs from the input (checks for duplicates)
    creates a metadata dataframe to keep track which IDs belong to which marker
    """
    output_path = ensure_path(output, "combined_db.faa", force = force)

    rows = []

    with open(output_path, "w") as out:
        for marker, fasta_path in sorted(db_fastas.items()):
            for header, sequence in stream_fasta(fasta_path):
                db_id = clean_fasta_header(header)

                out.write(f">{db_id}\n")
                out.write(f"{sequence}\n")

                rows.append({
                    "db_ID": db_id,
                    "marker": marker,
                    "db_fasta": fasta_path
                })

    db_metadata = pd.DataFrame(rows)

    duplicate_ids = db_metadata[db_metadata["db_ID"].duplicated(keep=False)]

    if not duplicate_ids.empty:
        raise ValueError(f"Duplicate IDs found in database: {duplicate_ids}")

    logger.info(f"Concatenated {len(db_fastas)} database FASTA file(s)")

    return output_path, db_metadata

def classification_output_new(
        query_path: str,
        bsr_path: str,
        yaml_files: dict[str, str],
        db_metadata: pd.DataFrame,
        output: str,
        force: bool = False):
    """
    classifies hits based off of marker-specific YAML cutoffs
    creates one output file that includes hits and non-hits
    for every protein and marker:
        check if there is a bsr hit for that marker
        classify using correct cutoffs
        write combined output
    """

    cutoffs_by_marker = {
        marker: load_yaml_cutoffs(yaml_path) for marker, yaml_path in yaml_files.items()
    }

    query_df = create_output_rows(query_path)

    marker_df = pd.DataFrame({
        "marker": sorted(yaml_files.keys())
    })

    result_df = query_df.merge(marker_df, how="cross")

    bsr_df = pd.read_csv(bsr_path, sep="\t")

    required_bsr_columns = [
        "qseqid",
        "sseqid",
        "pident",
        "qlen",
        "score",
        "max_score",
        "BSR"
    ]

    if not bsr_df.empty:
        for column in required_bsr_columns:
            if column not in bsr_df.columns:
                raise ValueError(f"Required column {column} not found in BSR table")

        bsr_df = bsr_df.rename(columns={
            "qseqid": "prot_ID",
            "sseqid": "db_ID",
        })

        db_marker_df = db_metadata[["db_ID", "marker"]].drop_duplicates()

        bsr_df = bsr_df.merge(db_marker_df, on="db_ID", how="left")

        missing_marker = bsr_df[bsr_df["marker"].isna()]

        if not missing_marker.empty:
            raise ValueError(
                "Some BSR hits could not be assigned to a marker"
            )

        bsr_df["classification"] = bsr_df.apply(
            lambda row: classify_by_yaml_cutoffs(
                max_score=row["max_score"],
                hit_score=row["score"],
                bsr_cutoff=cutoffs_by_marker[row["marker"]]["bsr_cutoff"],
                max_score_lower=cutoffs_by_marker[row["marker"]]["max_score_lower"],
                max_score_upper=cutoffs_by_marker[row["marker"]]["max_score_upper"]
            ),
            axis=1
        )

        bsr_df = bsr_df.sort_values(
            by=["prot_ID", "marker", "score"],
            ascending=[True, True, False])

        bsr_df = bsr_df.drop_duplicates(
            subset=["prot_ID", "marker"],
            keep="first"
        )

        result_df = result_df.merge(
            bsr_df,
            on=["prot_ID", "marker"],
            how="left"
        )

    result_df["classification"] = result_df["classification"].fillna("no_hit")

    """
    result_df["annotation"] = result_df.apply(
        lambda row:
            cutoffs_by_marker[row["marker"]]["protein_name"]
            if row["classification"] not in ["below_cutoff", "no_hit"]
            else pd.NA,
        axis=1
    )
    # if we want to keep the "annotation" columns: add it to first_columns
    """

    result_df = result_df.replace("no_hit", pd.NA)

    first_columns = [
        "prot_ID",
        "sequence_length",
        "marker",
        "classification",
    ]

    additional_columns = []

    for column in result_df.columns:
        if column not in first_columns:
            additional_columns.append(column)

    result_df = result_df[first_columns + additional_columns]

    output_path = ensure_path(output, "annotate.tsv", force=force)

    result_df.to_csv(output_path, sep="\t", index=False)

    logger.info(f"Classification complete! Hits saved to: {output_path}")

    return output_path

def annotate_plot_by_marker(
        bsr_path: str,
        yaml_files: dict[str, str],
        db_metadata: pd.DataFrame,
        output: str,
        svg: bool = False,
        force: bool = False):
    """
    creates a dictionary for all plots
    creates one plot per marker
    """

    plot_dir = Path(output) / "annotate_plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    bsr_df = pd.read_csv(bsr_path, sep="\t")

    if bsr_df.empty:
        logger.warning("No hits found in BSR table. Skipping all annotate plots.")
        return []

    required_bsr_columns = [
        "qseqid",
        "sseqid",
        "pident",
        "qlen",
        "score",
        "max_score",
        "BSR"
    ]

    for column in required_bsr_columns:
        if column not in bsr_df.columns:
            raise ValueError(f"Required column {column} not found in BSR table")

    bsr_df = bsr_df.rename(columns={
        "qseqid": "prot_ID",
        "sseqid": "db_ID"
    })

    db_marker_df = db_metadata[["db_ID", "marker"]].drop_duplicates()

    bsr_df = bsr_df.merge(db_marker_df, on="db_ID", how="left")

    missing_marker = bsr_df[bsr_df["marker"].isna()]

    if not missing_marker.empty:
        raise ValueError("Some BSR hits could not be assigned to a marker")

    plot_paths = []

    for marker, yaml_path in yaml_files.items():
        marker_bsr_df = bsr_df[bsr_df["marker"] == marker].copy()

        plot_path = annotate_plot_helper(
            bsr_df=marker_bsr_df,
            yaml_path=yaml_path,
            output=str(plot_dir),
            marker=marker,
            svg=svg,
            force=force
        )

        if plot_path is not None:
            plot_paths.append(plot_path)

    logger.info(f"Created {len(plot_paths)} annotate plot(s)")

    return plot_paths

def annotate(
        db_path: str,
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

    # collect marker-specific files
    db_fastas = collect_db_fastas(db_path)
    yaml_files = collect_yaml_files(db_path)

    check_db_yaml_matching(
        db_markers=set(db_fastas.keys()),
        yaml_markers=set(yaml_files.keys())
    )

    # concatenate DB FASTAS
    combined_db_path, db_metadata = concatenate_db_input(
        db_fastas=db_fastas,
        output=output,
        force=force,
    )

    # build the DIAMOND DB from combined DB FASTA
    diamond_db = build(
        seed_fasta=combined_db_path,
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
    annotate_tsv = classification_output_new(
        query_path=query_path,
        bsr_path=bsr_path,
        yaml_files=yaml_files,
        db_metadata=db_metadata,
        output=output,
        force=force
    )

    annotate_plot_paths = annotate_plot_by_marker(
        bsr_path=bsr_path,
        yaml_files=yaml_files,
        db_metadata=db_metadata,
        output=output,
        force=force
    )

    logger.info("Annotate workflow completed")

    return annotate_tsv