#!/usr/bin/env python3

from .util import ensure_path, read_fasta_to_dict, parse_protein_identifier, extract_cog_info

import pandas as pd
import logging
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable, get_cmap
from typing import Optional
from pathlib import Path
import yaml
import sys
import tarfile
import gzip

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.DEBUG,  # or logging.INFO if you want less verbosity
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    stream=sys.stdout
)
logging.getLogger('matplotlib').setLevel(logging.WARNING)
logging.getLogger("PIL").setLevel(logging.INFO)

# ===========================
# CUGO GFF PARSER
# ===========================
def get_tmhmm_data_for_file(tmhmm_tar: tarfile.TarFile, file_id: str) -> dict:
    tmhmm_filename = f"{file_id}_tmhmm_clean.gz"

    for member in tmhmm_tar.getmembers():
        if member.name.endswith(tmhmm_filename):
            file_obj = tmhmm_tar.extractfile(member)
            if file_obj is None:
                return {}

            try:
                content = gzip.decompress(file_obj.read()).decode('utf-8')
                file_tmhmm = {}

                for line in content.strip().split('\n'):
                    if line.startswith('prot_ID'):
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        prot_id = parts[0]
                        no_tmh = int(parts[1])
                        file_tmhmm[prot_id] = no_tmh

                return file_tmhmm
            except Exception:
                return {}

    return {}


def extract_gene_info_with_tmhmm(line: list, tmhmm_dict: dict) -> tuple:
    annotation = line[8].split(';')
    seqID = annotation[0].split('=')[1].replace('__', '_')
    COG_ID = annotation[1].split('=')[1] if len(annotation) > 1 else 'NA'
    parent = line[0]
    direction = line[6]
    gene_start = line[3]
    gene_end = line[4]
    nuc_length = abs(int(gene_end) - int(gene_start)) + 1
    aa_length = int(nuc_length / 3)

    no_tmh = tmhmm_dict.get(COG_ID, 0)

    return seqID, COG_ID, parent, direction, gene_start, gene_end, nuc_length, aa_length, no_tmh

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

def parse_single_gff_with_tmhmm(gff_content: str, tmhmm_dict: dict) -> list:
    reformat_data = []
    cugo_size = {}
    cugo_count = prev_direction = prev_parent = prev_feat_type = prev_cugo = 0
    cugo_size_count = 0
    prev_line = None
    cds_count = 0

    lines = gff_content.strip().split('\n')

    for line_count, line in enumerate(lines, 1):
        clean_line = line.strip().split("\t")

        if len(clean_line) != 9:
            continue

        feat_type = clean_line[2]
        if feat_type == "CDS":
            cds_count += 1

        if feat_type != "CDS" and prev_feat_type != "CDS":
            continue
        prev_feat_type = feat_type

        if prev_line is None:
            prev_line = clean_line
            continue

        seqID, COG_ID, parent, direction, gene_start, gene_end, nuc_length, aa_length, no_tmh = extract_gene_info_with_tmhmm(
            prev_line, tmhmm_dict)

        next_parent = clean_line[0]
        next_direction = clean_line[6]

        cugo, cugo_start, cugo_end, cugo_count, cugo_size_count = cugo_boundaries(
            direction, prev_direction, next_direction,
            parent, prev_parent, next_parent,
            cugo_count, cugo_size, cugo_size_count, prev_cugo
        )

        reformat_data.append([seqID, parent, gene_start, gene_end, nuc_length, aa_length,
                              direction, COG_ID, cugo, cugo_start, cugo_end, no_tmh])

        prev_line = clean_line
        prev_direction = direction
        prev_parent = parent
        prev_cugo = cugo

    if len(clean_line) == 9 and clean_line[2] == "CDS":
        seqID, COG_ID, parent, direction, gene_start, gene_end, nuc_length, aa_length, no_tmh = extract_gene_info_with_tmhmm(
            clean_line, tmhmm_dict)

        cugo, cugo_start, cugo_end, cugo_count, cugo_size_count = cugo_boundaries(
            direction, prev_direction, 0,
            parent, prev_parent, "",
            cugo_count, cugo_size, cugo_size_count, prev_cugo
        )

        reformat_data.append([seqID, parent, gene_start, gene_end, nuc_length, aa_length,
                              direction, COG_ID, cugo, cugo_start, cugo_end, no_tmh])

        cugo_size[cugo] = cugo_size_count


    final_data = []
    for row in reformat_data:
        seqID = row[0]
        parent_ID = row[1]
        aa_length = row[5]
        strand = row[6]
        COG_ID = row[7]
        cugo_number = row[8]
        no_tmh = row[11]

        final_data.append([seqID, parent_ID, aa_length, strand, COG_ID, cugo_number, no_tmh])

    return final_data


def get_file_id_from_gff_name(gff_name: str) -> str:
    base_name = Path(gff_name).name
    if '_cog.gff' in base_name:
        return base_name.split('_cog.gff')[0]
    return base_name


def parse(tar_gz_path: str, tmhmm_tar_path: str = None, output_dir: str = None, globdb_version: int = None,
          force: bool = False):
    output_path = ensure_path(output_dir, f"globdb_r{globdb_version}_cugo")

    file_count = 0
    columns = ["seqID", "parent_ID", "aa_length", "strand", "COG_ID", "CUGO_number", "no_TMH"]

    with open(output_path, 'w') as output_file:
        output_file.write('\t'.join(columns) + '\n')

        try:
            with tarfile.open(tar_gz_path, 'r:gz') as gff_tar:
                # Open TMHMM tar only if path is provided
                tmhmm_tar = None
                if tmhmm_tar_path:
                    tmhmm_tar = tarfile.open(tmhmm_tar_path, 'r:gz')

                try:
                    members = gff_tar.getmembers()
                    total_files = len([m for m in members if '_cog.gff' in m.name])

                    logger.info(f"Found {total_files} GFF files to process")
                    if tmhmm_tar_path:
                        logger.info("TMHMM data will be included")
                    else:
                        logger.info("TMHMM data will be skipped (no TMHMM tar file provided)")

                    for member in members:
                        if '_cog.gff' not in member.name:
                            continue

                        file_count += 1

                        file_id = get_file_id_from_gff_name(member.name)
                        file_tmhmm = {}

                        if tmhmm_tar:
                            tmhmm_data = get_tmhmm_data_for_file(tmhmm_tar, file_id)
                            if tmhmm_data:
                                file_tmhmm = tmhmm_data
                            else:
                                logger.warning(f"No TMHMM data found for {file_id}")

                        file_obj = gff_tar.extractfile(member)
                        if file_obj is None:
                            logger.warning(f"Could not extract {member.name}")
                            continue

                        try:
                            if member.name.endswith('.gz'):
                                content = gzip.decompress(file_obj.read()).decode('utf-8')
                            else:
                                content = file_obj.read().decode('utf-8')

                            file_data = parse_single_gff_with_tmhmm(content, file_tmhmm)

                            for row in file_data:
                                output_file.write('\t'.join(map(str, row)) + '\n')

                        except Exception as e:
                            logger.error(f"Error processing {member.name}: {str(e)}")
                            continue

                        if file_count % 1000 == 0:
                            logger.info(f"Processed {file_count}/{total_files} files")

                finally:
                    if tmhmm_tar:
                        tmhmm_tar.close()

        except Exception as e:
            logger.error(f"Error processing tar.gz file: {str(e)}")
            raise

    logger.info(f"Successfully processed {file_count} files. Output saved to {output_path}")

# ===========================
# CUGO CONTEXT PARSER
# ===========================

def context(protein_ids: Optional[str],
            cugo_path: Path,
            cugo_range: int,
            output_dir: str,
            protein_name: str,
            force: bool = False,
            fasta_path: Optional[str] = None):
    log_file = ensure_path(output_dir, f"{protein_name}_missing_files.log", force=force)
    logging.basicConfig(
        filename=log_file,
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    # Load protein identifiers
    if fasta_path:
        sequences = read_fasta_to_dict(fasta_path)
        protein_identifiers = list(sequences.keys())
    elif protein_ids:
        with open(protein_ids, 'r') as id_file:
            protein_identifiers = [line.strip() for line in id_file if line.strip()]
    else:
        logger.error("Either 'fasta_path' or 'protein_ids' must be provided.")
        raise ValueError("You must provide either a FASTA file or a list of protein IDs.")

    output_file = ensure_path(output_dir, f"{protein_name}_context.tsv", force=force)

    with gzip.open(cugo_path, "rt") as f:
        df = pd.read_csv(f, sep="\t", na_filter=False)

    # Set up results container
    results = None
    protein_set = set(protein_identifiers)

    # Filter for proteins of interest
    target_df = df[df["seqID"].isin(protein_set)]

    if target_df.empty:
        logging.warning("No target proteins found in the CUGO file.")
        return

    for _, row in target_df.iterrows():
        id = row["seqID"]
        target_cugo = int(row["CUGO_number"])
        target_parent = row["parent_ID"]
        target_strand = row["strand"]

        # Get context window for this protein
        parent_df = df[df["parent_ID"] == target_parent]
        cugo_context = parent_df[
            (parent_df["CUGO_number"] >= (target_cugo - cugo_range)) &
            (parent_df["CUGO_number"] <= (target_cugo + cugo_range))
        ]
        if target_strand == "-":
            cugo_context = cugo_context.iloc[::-1]

        cugo_context = cugo_context.reset_index(drop=True)

        # Find the index of the target
        target_index = cugo_context.index[cugo_context.seqID == id].item()
        cugo_context.index = cugo_context.index - target_index
        cugo_context = cugo_context.transpose()

        # Store
        if results is None:
            results = cugo_context
        else:
            results = pd.concat([results, cugo_context], join="outer", axis=0)

    # Save results
    if results is not None:
        results = results.reset_index().rename(columns={"index": "feat_type"})
        results.to_csv(output_file, sep="\t", index=False)


# ======================================
# FUNCTION DEFINITIONS FOR CUGO PLOTTING
# ======================================
def load_cugo_context(path: str):
    return pd.read_csv(path, sep='\t', na_values='', keep_default_na=False)

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
        title: str = "Top COGs per Position",
        show: bool = False,
        ax: plt.Axes = None
):
    cont = load_cugo_context(cugo_path)
    extract = extract_cog_info(cont, "COG_ID")

    script_dir = Path(__file__).resolve().parent
    color_yaml_path = script_dir / "cog_colors.yaml"

    with open(color_yaml_path, "r") as f:
        cog_color_map = yaml.safe_load(f)

    flank_cols = [str(i) for i in range(flank_lower, flank_upper + 1)]
    cugo_df = extract[flank_cols]
    positions = [int(col) for col in cugo_df.columns]

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

    subtick_width = 0.8 / top_n
    subtick_offset = (top_n - 1) * subtick_width / 2

    x_pos, y_values, cog_labels, point_colors = [], [], [], []
    all_xticks, all_xlabels = [], []

    if ax is None:
        figsize = (max(8, len(positions) * 0.8), 7)
        fig, ax = plt.subplots(figsize=figsize)

    for pos_idx, pos in enumerate(positions):
        for rank in range(top_n):
            cog_id = id_df.loc[rank, str(pos)]
            count = count_df.loc[rank, str(pos)]

            x = pos + (rank * subtick_width) - subtick_offset
            x_pos.append(x)
            y_values.append(count)
            all_xticks.append(x)

            if pd.isna(cog_id):
                cog_labels.append("nan")
                point_colors.append("#ffffff")
                all_xlabels.append("nan")
            else:
                cog_labels.append(cog_id)
                point_colors.append(cog_color_map.get(cog_id, "#aaaaaa"))
                all_xlabels.append(cog_id)

    ax.scatter(
        x_pos,
        y_values,
        color=point_colors,
        s=150,
        zorder=3,
        edgecolors='black',
        linewidths=0.7
    )

    ax.set_xticks(all_xticks)
    ax.set_xticklabels(all_xlabels, rotation=90, fontsize=14, ha='center')

    for pos in positions[:-1]:
        boundary = pos + 0.5
        ax.axvline(boundary, color='gray', linestyle='--', linewidth=0.7, zorder=1)

    max_y = max([y for y in y_values if not pd.isna(y)], default=0)
    label_offset = max_y * 0.05

    for pos in positions:
        ax.text(pos, max_y + label_offset, str(pos), ha='center', va='bottom',
                fontsize=16, fontweight='bold')

    ax.set_ylim(-max_y * 0.05, max_y * 1.15)
    ax.set_xlim(positions[0] - 0.5, positions[-1] + 0.5)

    ax.tick_params(axis='y', labelsize=14)
    ax.set_ylabel("Count", fontsize=16)
    ax.set_title(title, fontsize=18, fontweight='bold')

    if show:
        plt.tight_layout()
        plt.show()

    return positions, [pos + 0.5 for pos in positions[:-1]]

def plot_size_per_position(cugo_path: str,
                           flank_lower: int,
                           flank_upper: int,
                           title: str = "aa_length Density per Position",
                           show: bool = False,
                           ax: Optional[plt.Axes] = None,
                           pos_boundaries: Optional[list] = None,
                           bin_width: int = 10,
                           y_range: int = None
                           ):
    cont = pd.read_csv(cugo_path, sep='\t', na_values='', keep_default_na=False)
    extract = cont.loc[cont["feat_type"] == "aa_length"]

    flank_cols = [str(i) for i in range(flank_lower, flank_upper + 1)]
    cugo_df = extract[flank_cols]
    positions = [int(col) for col in cugo_df.columns]

    all_lengths = cugo_df.values.flatten()
    all_lengths = all_lengths[~pd.isna(all_lengths)].astype(float)
    max_len = all_lengths.max()
    bin_edges = np.arange(0, max_len + bin_width, bin_width)
    n_bins = len(bin_edges) - 1

    heat_data = np.zeros((n_bins, len(flank_cols)))
    position_counts = []

    for col_idx, col in enumerate(flank_cols):
        values = cugo_df[col].dropna().astype(float)
        hist, _ = np.histogram(values, bins=bin_edges)
        heat_data[:, col_idx] = hist
        position_counts.append(len(values))

    if ax is None:
        fig, ax = plt.subplots()

    cmap = get_cmap("Blues")
    norm = Normalize(vmin=0, vmax=heat_data.max())

    rect_width = 0.8
    for col_idx, pos in enumerate(positions):
        for y in range(heat_data.shape[0]):
            val = heat_data[y, col_idx]
            color = cmap(norm(val))
            ax.add_patch(plt.Rectangle(
                (pos - rect_width / 2, y), rect_width, 1, color=color, linewidth=0
            ))

    if pos_boundaries is not None:
        for boundary in pos_boundaries:
            ax.axvline(boundary, color='gray', linestyle='--', linewidth=0.7, zorder=1)

    ax.set_xlim(positions[0] - 0.5, positions[-1] + 0.5)
    xtick_labels = [f"{pos}\n(n={count})" for pos, count in zip(positions, position_counts)]
    ax.set_xticks(positions)
    ax.set_xticklabels(xtick_labels, fontsize=14)

    if y_range:
        n_bins = int(y_range / bin_width)

    ax.set_ylim(0, n_bins)

    tick_indices = []
    tick_labels = []

    for i in range(n_bins):
        bin_end = int(bin_edges[i + 1])
        if bin_end % 100 == 0:
            tick_indices.append(i)
            tick_labels.append(f"{bin_end}")

    if not tick_indices:
        tick_step = max(1, n_bins // 10)
        tick_indices = list(range(0, n_bins, tick_step))
        tick_labels = [f"{int(bin_edges[i + 1])}" for i in tick_indices]

    ax.set_yticks(np.array(tick_indices) + 0.5)
    ax.set_yticklabels(tick_labels, fontsize=16)

    ax.set_xlabel("Position", fontsize=18)
    ax.set_ylabel(f"length (bin size: {bin_width})", fontsize=18)
    ax.set_title(title, fontsize=20)

    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])

    plt.tight_layout()
    if show:
        plt.show()



def cugo_plot(cugo_path: str,
              flank_lower: int,
              flank_upper: int,
              top_n: int,
              cugo: bool = False,
              size: bool = False,
              all_plots: bool = False,
              bin_width: int = 10,
              y_range: int = None
              ):
    if cugo:
        if not top_n:
            logger.error("top_n is required if cugo is True")
        else:
            plot_top_cogs_per_position(cugo_path, flank_lower, flank_upper, top_n, show=True)
    if size:
        plot_size_per_position(cugo_path, flank_lower, flank_upper, show=True, bin_width=bin_width)
    if all_plots:
        if not top_n:
            logger.error("top_n is required if cugo is True")
        else:
            width = max(8, (flank_upper - flank_lower + 1) * top_n * 0.6)
            figsize = (width, 12)
            fig, (ax1, ax2) = plt.subplots(2, 1,
                                           figsize=figsize,
                                           gridspec_kw={'height_ratios': [3, 2]})

            pos_centers, pos_boundaries = plot_top_cogs_per_position(
                cugo_path=cugo_path,
                flank_lower=flank_lower,
                flank_upper=flank_upper,
                top_n=top_n,
                title="Top COGs per Position",
                ax=ax1,
                show=False
            )

            plot_size_per_position(
                cugo_path=cugo_path,
                flank_lower=flank_lower,
                flank_upper=flank_upper,
                title="aa_length Density per Position",
                ax=ax2,
                show=False,
                pos_boundaries=pos_boundaries,
                bin_width=bin_width,
                y_range=y_range
            )

            ax1.tick_params(labelbottom=True)


            plt.tight_layout()
            plt.subplots_adjust(hspace=0.4)
            plt.show()