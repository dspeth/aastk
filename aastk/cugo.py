#!/usr/bin/env python3

from .util import ensure_path

import pandas as pd
import logging
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import yaml

# ========================================
# FUNCTION DEFINITIONS FOR CUGO GFF PARSER
# ========================================

# --- Gene Info Extraction ---
def extract_gene_info(line):
    """
    Extracts gene info from a GFF line
    Args:
        line (list): A tab-separated line from GFF file split into fields

    Returns:
        tuple: (seqID, COG_ID, parent, direction, gene_start, gene_end, nuc_length, aa_length)
        seqID (str): Sequence identifier, cleaned by replacing ___ with _
        COG_ID (str): COG identifier or "NA" if not present
        parent (str): Parent sequence/contig identifier
        direction (str): Strand direction ("+" or "-")
        gene_start (str): Start position of the gene
        gene_end (str): End position of the gene
        nuc_length (int): Nucleotide length of the gene
        aa_length (int): Amino acid length (nucleotide length / 3)
    """
    annotation = line[8].split(';')
    seqID = annotation[0].split('=')[1].replace('__', '_')
    COG_ID = annotation[1].split('=')[1] if len(annotation) > 1 else 'NA'
    parent = line[0]
    direction = line[6]
    gene_start = line[3]
    gene_end = line[4]
    nuc_length = abs(int(gene_end) - int(gene_start)) + 1
    aa_length = int(nuc_length / 3)
    return seqID, COG_ID, parent, direction, gene_start, gene_end, nuc_length, aa_length

def cugo_boundaries(direction: str,
                    prev_direction: str,
                    next_direction: str,
                    parent: str,
                    prev_parent: str,
                    next_parent: str,
                    cugo_count: int,
                    cugo_size: dict,
                    cugo_size_count: int,
                    prev_cugo: int = None):
    cugo_start, cugo_end = "NA", "NA"

    # middle of CUGO
    if (direction == prev_direction == next_direction) and (parent == prev_parent == next_parent):
        cugo = prev_cugo
        cugo_size_count += 1

    # start of contig/scaffold
    elif parent != prev_parent:
        cugo = cugo_count
        cugo_size_count = 1

        if direction == "+":
            cugo_start = "sequence_edge"
            if parent != next_parent:
                cugo_end = "sequence_edge"
                cugo_size[cugo] = cugo_size_count
                cugo_size_count = 0
                cugo_count += 1
            elif direction != next_direction:
                cugo_end = "strand_change"
                cugo_size[cugo] = cugo_size_count
                cugo_size_count = 0
                cugo_count += 1
        else:  # direction == "-"
            cugo_end = "sequence_edge"
            if parent != next_parent:
                cugo_start = "sequence_edge"
                cugo_size[cugo] = cugo_size_count
                cugo_size_count = 0
                cugo_count += 1
            elif direction != next_direction:
                cugo_start = "strand_change"
                cugo_size[cugo] = cugo_size_count
                cugo_size_count = 0
                cugo_count += 1

    # end of contig/scaffold
    elif parent != next_parent:
        if direction == prev_direction:
            cugo = prev_cugo
            cugo_size_count += 1
            if direction == "+":
                cugo_end = "sequence_edge"
            else:
                cugo_start = "sequence_edge"
            cugo_size[cugo] = cugo_size_count
            cugo_size_count = 0
            cugo_count += 1
        else:
            cugo = cugo_count
            cugo_size_count = 1
            if direction == "+":
                cugo_start = "strand_change"
                cugo_end = "sequence_edge"
            else:
                cugo_start = "sequence_edge"
                cugo_end = "strand_change"
            cugo_size[cugo] = cugo_size_count
            cugo_size_count = 0
            cugo_count += 1

    # Strand change
    elif direction != prev_direction or direction != next_direction:
        if direction != prev_direction and direction == next_direction:
            cugo = cugo_count
            cugo_size_count = 1
            if direction == "+":
                cugo_start = "strand_change"
            else:
                cugo_end = "strand_change"
        elif direction == prev_direction and direction != next_direction:
            cugo = prev_cugo
            cugo_size_count += 1
            if direction == "+":
                cugo_end = "strand_change"
            else:
                cugo_start = "strand_change"
            cugo_size[cugo] = cugo_size_count
            cugo_size_count = 0
            cugo_count += 1
        else:  # direction != prev_direction and direction != next_direction
            cugo = cugo_count
            cugo_size_count = 1
            cugo_start = cugo_end = "strand_change"
            cugo_size[cugo] = cugo_size_count
            cugo_size_count = 0
            cugo_count += 1

    return cugo, cugo_start, cugo_end, cugo_count, cugo_size_count

def process_last_gene(line,
                      prev_direction: str,
                      prev_parent: str,
                      cugo_count: int,
                      cugo_size_count: int,
                      prev_cugo: int,
                      cugo_size: dict):
    seqID, COG_ID, parent, direction, gene_start, gene_end, nuc_length, aa_length = extract_gene_info(line)
    cugo_start, cugo_end = "NA", "NA"

    if (direction == prev_direction) and (parent == prev_parent):
        cugo = prev_cugo
        cugo_size_count += 1
        cugo_size[cugo] = cugo_size_count
        if direction == "+":
            cugo_end = "sequence_edge"
        else:
            cugo_start = "sequence_edge"
    elif parent != prev_parent:
        cugo = cugo_count
        cugo_size_count = 1
        cugo_size[cugo] = cugo_size_count
        cugo_start = cugo_end = "sequence_edge"
    elif direction != prev_direction:
        cugo = cugo_count
        cugo_size_count = 1
        cugo_size[cugo] = cugo_size_count
        if direction == "+":
            cugo_start, cugo_end = "strand_change", "sequence_edge"
        else:
            cugo_start, cugo_end = "sequence_edge", "strand_change"

    return [seqID, parent, gene_start, gene_end, nuc_length, aa_length,
            direction, COG_ID, cugo, cugo_start, cugo_end], cugo_size

# ===========================
# CUGO GFF PARSER
# ===========================
def parse(gff_file_path: str,
          output_dir: str):
    reformat_data = []
    cugo_size = {}
    cugo_count = prev_direction = prev_parent = prev_feat_type = prev_cugo = 0
    cugo_size_count = 0
    prev_line = None
    cds_count = 0

    # also implement reading directly from gzip - check if unpacked file is there, otherwise read
    # from gzip directly;
    # name eventual output file
    identifier = gff_file_path.split('.')[0].rsplit("_", 1)[0] + '_cugo.tsv'
    output_path = ensure_path(output_dir, identifier)

    try:
        with open(gff_file_path, "r") as GFF:
            for line_count, line in enumerate(GFF, 1):
                clean_line = line.strip().split("\t")

                # skip invalid lines and non-CDS features
                if len(clean_line) != 9:
                    continue

                feat_type = clean_line[2]
                if feat_type == "CDS":
                    cds_count += 1

                if feat_type != "CDS" and prev_feat_type != "CDS":
                    continue
                prev_feat_type = feat_type

                # initialize with first valid line
                if prev_line is None:
                    prev_line = clean_line
                    continue

                # extract information from current line
                seqID, COG_ID, parent, direction, gene_start, gene_end, nuc_length, aa_length = extract_gene_info(
                    prev_line)

                # get context
                next_parent = clean_line[0]
                next_direction = clean_line[6]

                # determine CUGO membership
                cugo, cugo_start, cugo_end, cugo_count, cugo_size_count = cugo_boundaries(
                    direction, prev_direction, next_direction,
                    parent, prev_parent, next_parent,
                    cugo_count, cugo_size, cugo_size_count, prev_cugo
                )

                # add to results
                reformat_data.append([seqID, parent, gene_start, gene_end, nuc_length, aa_length,
                                      direction, COG_ID, cugo, cugo_start, cugo_end])

                # update for next iteration
                prev_line = clean_line
                prev_direction = direction
                prev_parent = parent
                prev_cugo = cugo

            # process last line if it's valid
            if len(clean_line) == 9 and clean_line[2] == "CDS":
                last_line, cugo_size = process_last_gene(
                    clean_line, prev_direction, prev_parent,
                    cugo_count, cugo_size_count, prev_cugo, cugo_size
                )
                reformat_data.append(last_line)

        # create dataframe
        gff_cugo = pd.DataFrame(
            reformat_data,
            columns=["seqID", "parent_ID", "gene_start", "gene_end", "nuc_length",
                     "aa_length", "strand", "COG_ID", "CUGO_number", "CUGO_start", "CUGO_end"]
        )
        gff_cugo["CUGO_size"] = gff_cugo["CUGO_number"].map(cugo_size)

        gff_cugo.to_csv(output_path, sep="\t", index=False)

        return gff_cugo

    except Exception as e:
        raise Exception(f"Error processing GFF file: {str(e)}")

# ===========================
# CUGO CONTEXT PARSER
# ===========================
# combine with plotting - no intermediate files needed
def context(protein_ids: str, cugo_dir: str, tmhmm_dir: str, cugo_range: int, output_dir: str, dataset: str):
    log_file = ensure_path(output_dir, f"{dataset}_missing_files.log")
    logging.basicConfig(
        filename=log_file,
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    with open(protein_ids, 'r') as id_file:
        protein_identifiers = [line.strip() for line in id_file if line.strip()]
        clean_identifiers = [protein_identifier.rsplit("___", 1)[0] for protein_identifier in protein_identifiers]

    output_file = ensure_path(output_dir, dataset + "_context.tsv")

    missing_files = []
    results = None

    for index, id in enumerate(protein_identifiers):
        cugo_tab_path = Path((cugo_dir) / f"{clean_identifiers[index]}_cugo.gz")

        # Check if CUGO file exists
        if not cugo_tab_path.exists():
            logging.warning(f"Missing CUGO file: {cugo_tab_path}")
            missing_files.append(str(cugo_tab_path))
            continue

        # Try to read the gzipped file
        try:
            genome_cugo_df = pd.read_csv(cugo_tab_path, compression='gzip', sep='\t', na_filter=False)
        except Exception as e:
            logging.error(f"Error reading CUGO file {cugo_tab_path}: {str(e)}")
            missing_files.append(str(cugo_tab_path))
            continue

        # Process TMHMM file if provided
        if tmhmm_dir:
            tmhmm_file = Path(tmhmm_dir) / f"{clean_identifiers[index]}_tmhmm_clean"
            if tmhmm_file.exists():
                genome_tmhmm_df = pd.read_csv(tmhmm_file, sep="\t", na_filter=False)
                parse_df = pd.merge(genome_cugo_df, genome_tmhmm_df, on="prot_ID", how="left")
            else:
                parse_df = genome_cugo_df
        else:
            parse_df = genome_cugo_df

        # in the CUGO tab we only care about our protein
        target_select = parse_df[parse_df['seqID'] == id]
        target_cugo = target_select["CUGO_number"].iloc[0]
        target_parent = target_select["parent_ID"].iloc[0]
        target_strand = target_select["strand"].iloc[0]

        parent_df = parse_df[parse_df["parent_ID"] == target_parent]
        cugo_context = parent_df[(parent_df["CUGO_number"] >= (target_cugo - cugo_range)) &
                                 (parent_df["CUGO_number"] <= (target_cugo + cugo_range))]
        if target_strand == "-":
            cugo_context = cugo_context.iloc[::-1]
        cugo_context = cugo_context.reset_index(drop=True)

        # Find the index of our target protein
        target_index = cugo_context.index[cugo_context.seqID == id].item()
        cugo_context.index = cugo_context.index - target_index
        # Transpose the data
        cugo_context = cugo_context.transpose()

        # Store this protein's context data
        if results is None:
            # First protein - initialize results
            results = cugo_context
        else:
            # Concatenate with previous results
            # Use outer join to handle all possible columns from different proteins
            results = pd.concat([results, cugo_context], join="outer", axis=0)

        # If we have any results, save
    if results is not None:
        # Reset index and rename the index column
        results = results.reset_index().rename(columns={"index": "feat_type"})

        # Save the results
        results.to_csv(output_file, sep="\t", index=False)

# ======================================
# FUNCTION DEFINITIONS FOR CUGO PLOTTING
# ======================================
def load_cugo_context(path: str):
    return pd.read_csv(path, sep='\t', na_values='', keep_default_na=False)

def extract_cog_info(df: pd.DataFrame, feat_type=str):
    return df.loc[df["feat_type"] == feat_type]

def extract_flanking_window(df: pd.DataFrame, lower: int, upper: int):
    flank_list = list(range(lower, upper + 1))
    flank_str = [str(i) for i in flank_list]
    return df[flank_str]

def top_context(df: pd.DataFrame, top_n: int):
    top_ids = [[] for _ in range(top_n)]
    top_counts = [[] for _ in range(top_n)]

    for col in df.columns:
        counts = df[col].value_counts()
        for i in range(top_n):
            if i < len(counts):
                top_ids[i].append(counts.index[i])
                top_counts[i].append(counts.iloc[i])
            else:
                top_ids[i].append(np.nan)
                top_counts[i].append(np.nan)

    id_df = pd.DataFrame(top_ids, columns=df.columns)
    count_df = pd.DataFrame(top_counts, columns=df.columns)
    return id_df, count_df


def plot_top_cogs_per_position(
        cugo_path: str,
        flank_lower: int,
        flank_upper: int,
        top_n: int,
        title: str = "Top COGs per Position"
):
    """
    Plot the top N COGs per position from a CUGO context TSV file.

    Parameters:
    - cugo_path: Path to the CUGO context TSV file
    - flank_lower: Start of flanking window (inclusive)
    - flank_upper: End of flanking window (inclusive)
    - top_n: Number of top COGs to plot per position
    - title: Plot title string
    """

    # --- Load and prepare data ---
    cont = pd.read_csv(cugo_path, sep='\t', na_values='', keep_default_na=False)
    extract = cont.loc[cont["feat_type"] == "COG_ID"]

    script_dir = Path(__file__).resolve().parent
    color_yaml_path = script_dir / "cog_colors.yaml"

    with open(color_yaml_path, "r") as f:
        cog_color_map = yaml.safe_load(f)

    flank_cols = [str(i) for i in range(flank_lower, flank_upper + 1)]
    cugo_df = extract[flank_cols]

    # Calculate top N COGs per position
    top_ids = [[] for _ in range(top_n)]
    top_counts = [[] for _ in range(top_n)]

    for col in cugo_df.columns:
        counts = cugo_df[col].value_counts()
        for i in range(top_n):
            if i < len(counts):
                top_ids[i].append(counts.index[i])
                top_counts[i].append(counts.iloc[i])
            else:
                top_ids[i].append(np.nan)
                top_counts[i].append(np.nan)

    id_df = pd.DataFrame(top_ids, columns=cugo_df.columns)
    count_df = pd.DataFrame(top_counts, columns=cugo_df.columns)

    positions = id_df.columns
    plot_info = {
        pos: [
            {rank + 1: (id_df.loc[rank, pos], count_df.loc[rank, pos])}
            for rank in range(top_n)
        ]
        for pos in positions
    }

    point_spacing = 0.5 * top_n
    position_spacing = 0.5 * top_n

    x_pos, y_values, cog_labels, point_colors = [], [], [], []
    pos_boundaries, pos_centers = [], []

    current_x = 0
    for pos_idx, pos in enumerate(positions):
        pos_width = point_spacing * (top_n - 1)
        pos_center = current_x + pos_width / 2
        pos_centers.append(pos_center)

        if pos_idx == 0:
            pos_boundaries.append(current_x - point_spacing / 2)

        for rank in range(top_n):
            cog_id, count = list(plot_info[pos][rank].values())[0]
            x = current_x + rank * point_spacing
            x_pos.append(x)
            y_values.append(count)

            if pd.isna(cog_id):
                cog_labels.append("nan")
                point_colors.append("#ffffff")
            else:
                cog_labels.append(cog_id)
                point_colors.append(cog_color_map.get(cog_id, "#aaaaaa"))

        current_x += pos_width + position_spacing
        pos_boundaries.append(current_x - position_spacing / 2)

    figsize = (max(8, len(positions) * top_n * 0.6), 7)

    fig, ax = plt.subplots(figsize=figsize)

    ax.scatter(
        x_pos,
        y_values,
        color=point_colors,
        s=150,
        zorder=3,
        edgecolors='black',
        linewidths=0.7
    )

    ax.set_xticks(x_pos)
    ax.set_xticklabels(cog_labels, rotation=90, fontsize=14, ha='center')

    for boundary in pos_boundaries:
        ax.axvline(boundary, color='gray', linestyle='--', linewidth=0.7, zorder=1)

    max_y = max([y for y in y_values if not pd.isna(y)], default=0)
    label_offset = max_y * 0.05

    for center, pos in zip(pos_centers, positions):
        ax.text(center, max_y + label_offset, pos, ha='center', va='bottom',
                fontsize=16, fontweight='bold')

    ax.set_ylim(-max_y * 0.05, max_y * 1.15)
    ax.set_xlim(pos_boundaries[0] - point_spacing / 2, pos_boundaries[-1] + point_spacing / 2)

    ax.tick_params(axis='y', labelsize=14)
    ax.set_ylabel("Count", fontsize=16)
    ax.set_title(title, fontsize=18, fontweight='bold')

    plt.tight_layout()
    plt.subplots_adjust(left=0.05, bottom=0.25)

def cugo_plot(cugo_path: str,
        flank_lower: int,
        flank_upper: int,
        top_n: int):
    plot_top_cogs_per_position(cugo_path, flank_lower, flank_upper, top_n)
    plt.show()

### maybe plot transmembrane segments parts to distinguish the numbers contributing to specks;
### plotting should be on transposed context dataframe; THIS SHOULD BE SEPARATE FILES (length, cog, tmhmm) Columns by position, protein ID in column 1 and needs also genomeID column - maybe think about dynamic dataframes;
### processing slice by slice might be the quickest
### for - strand, transpose and reverse - depending on length, previous lines might need padding
### NA absence of annotation value (are genes, but don't have functional annotations), NaN no value


