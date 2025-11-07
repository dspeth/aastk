#!/usr/bin/env python3

from .util import *

import pandas as pd
import logging
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable, get_cmap
from typing import Optional
from pathlib import Path
import yaml
import gzip
import sqlite3
from tqdm import tqdm
import subprocess
from collections import defaultdict
import csv
from concurrent.futures import ThreadPoolExecutor, as_completed


logger = logging.getLogger(__name__)

# ===========================
# CUGO GFF PARSER
# ===========================
def setup_database(db_path: str) -> sqlite3.Connection:
    conn = sqlite3.connect(db_path)

    conn.execute('''
        CREATE TABLE IF NOT EXISTS all_data (
            seqID TEXT PRIMARY KEY,
            parent_ID TEXT,
            aa_length INTEGER,
            strand TEXT, 
            COG_ID TEXT,
            KEGG_ID TEXT,
            Pfam_ID TEXT,
            cugo_number INTEGER,
            no_tmh INTEGER
            )
        ''')

    conn.execute('CREATE INDEX IF NOT EXISTS idx_all_seqid ON all_data(seqID)')
    conn.execute('CREATE INDEX IF NOT EXISTS idx_parent_id ON all_data(parent_ID)')

    conn.commit()
    return conn


def extract_gene_info(line: list) -> tuple:
    """
    Extracts gene information from a GFF/GTF annotation line.

    Args:
        line (list): parsed annotation line with genomic feature data

    Returns:
        seqID (str): sequence identifier with underscores normalized
        annotation_ID (str): annotation identifier or 'NA' if not found
        parent (str): chromosome/contig name
        direction (str): strand direction (+/-)
        gene_start (str): start position
        gene_end (str): end position
        nuc_length (int): nucleotide sequence length
        aa_length (int): amino acid sequence length
    """
    # extract basic genomic coordinates and features
    parent = line[0]  # chromosome/contig name
    direction = line[6]  # strand direction (+/-)
    gene_start = line[3]  # start position
    gene_end = line[4]  # end position
    attributes = line[8]

    # parse attributes field
    attr_dict = {}
    for attr in attributes.split(';'):
        if '=' in attr:
            key, value = attr.split('=', 1)
            attr_dict[key] = value

    seqID = attr_dict.get('ID', '')
    annotation_ID = attr_dict.get('Name', '')


# calculate sequence lengths
    nuc_length = abs(int(gene_end) - int(gene_start)) + 1  # nucleotide length
    aa_length = int(nuc_length / 3)  # amino acid length
    return seqID, annotation_ID, parent, direction, gene_start, gene_end, nuc_length, aa_length


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
    """
    Determines CUGO (Co-oriented Unidirectional Gene Order) boundaries and classifications.

    Args:
        direction (str): current gene strand direction (+/-)
        prev_direction (str): previous gene strand direction (+/-)
        next_direction (str): next gene strand direction (+/-)
        parent (str): current gene's chromosome/contig name
        prev_parent (str): previous gene's chromosome/contig name
        next_parent (str): next gene's chromosome/contig name
        cugo_count (int): running count of CUGO regions
        cugo_size (dict): mapping of CUGO IDs to their sizes
        cugo_size_count (int): running count of genes in current CUGO
        prev_cugo (int, optional): previous CUGO identifier

    Returns:
        cugo (int): current CUGO identifier
        cugo_start (str): start boundary type ('sequence_edge', 'strand_change', or 'NA')
        cugo_end (str): end boundary type ('sequence_edge', 'strand_change', or 'NA')
        cugo_count (int): updated CUGO count
        cugo_size_count (int): updated gene count for current CUGO
    """
    cugo, cugo_start, cugo_end = 'NA', 'NA', 'NA'

    # middle of cugo - same direction and parent for all three genes
    if (direction == prev_direction == next_direction) and (parent == prev_parent == next_parent):
        cugo = prev_cugo
        cugo_size_count += 1

    # start of contig/scaffold - different parent from previous gene
    elif parent != prev_parent:
        cugo = cugo_count
        cugo_size_count = 1

        if direction == '+':
            cugo_start = 'sequence_edge'
            if parent != next_parent:
                # single gene contig
                cugo_end = 'sequence_edge'
                cugo_size[cugo] = cugo_size_count
                cugo_size_count = 0
                cugo_count += 1
            elif direction != next_direction:
                # strand change at start of contig
                cugo_end = 'strand_change'
                cugo_size[cugo] = cugo_size_count
                cugo_size_count = 0
                cugo_count += 1
        else:  # direction == '-'
            cugo_end = 'sequence_edge'
            if parent != next_parent:
                # single gene contig
                cugo_start = 'sequence_edge'
                cugo_size[cugo] = cugo_size_count
                cugo_size_count = 0
                cugo_count += 1
            elif direction != next_direction:
                # strand change at start of contig
                cugo_start = 'strand_change'
                cugo_size[cugo] = cugo_size_count
                cugo_size_count = 0
                cugo_count += 1

    # end of contig/scaffold - different parent from next gene
    elif parent != next_parent:
        if direction == prev_direction:
            # continuing previous cugo
            cugo = prev_cugo
            cugo_size_count += 1
            if direction == '+':
                cugo_end = 'sequence_edge'
            else:
                cugo_start = 'sequence_edge'
            cugo_size[cugo] = cugo_size_count
            cugo_size_count = 0
            cugo_count += 1
        else:
            # new cugo at end of contig
            cugo = cugo_count
            cugo_size_count = 1
            if direction == '+':
                cugo_start = 'strand_change'
                cugo_end = 'sequence_edge'
            else:
                cugo_start = 'sequence_edge'
                cugo_end = 'strand_change'
            cugo_size[cugo] = cugo_size_count
            cugo_size_count = 0
            cugo_count += 1

    # strand change - different direction from previous or next gene
    elif direction != prev_direction or direction != next_direction:
        if direction != prev_direction and direction == next_direction:
            # start of new cugo due to strand change
            cugo = cugo_count
            cugo_size_count = 1
            if direction == '+':
                cugo_start = 'strand_change'
            else:
                cugo_end = 'strand_change'
        elif direction == prev_direction and direction != next_direction:
            # end of current cugo due to strand change
            cugo = prev_cugo
            cugo_size_count += 1
            if direction == '+':
                cugo_end = 'strand_change'
            else:
                cugo_start = 'strand_change'
            cugo_size[cugo] = cugo_size_count
            cugo_size_count = 0
            cugo_count += 1
        else:  # direction != prev_direction and direction != next_direction
            # isolated gene with different direction from both neighbors
            cugo = cugo_count
            cugo_size_count = 1
            cugo_start = cugo_end = 'strand_change'
            cugo_size[cugo] = cugo_size_count
            cugo_size_count = 0
            cugo_count += 1

    return cugo, cugo_start, cugo_end, cugo_count, cugo_size_count

def parse_single_gff(gff_content: str) -> list:
    """
    Parses GFF content and extracts gene information with CUGO boundary analysis.

    Args:
        gff_content (str): GFF format file content as string

    Returns:
        list: List of processed gene records, each containing:
            - seqID (str): gene/sequence identifier
            - parent_ID (str): chromosome/contig identifier
            - aa_length (int): amino acid sequence length
            - strand (str): gene strand direction (+/-)
            - annotation ID (str): functional category identifier
            - cugo_number (int): CUGO region identifier
    """
    reformat_data = []
    cugo_size = {}
    cugo_count = prev_direction = prev_parent = prev_feat_type = prev_cugo = 0
    cugo_size_count = 0
    prev_line = None
    cds_count = 0

    # split content into individual lines for processing
    lines = gff_content.strip().split('\n')

    for line_count, line in enumerate(lines, 1):
        # parse tab-delimited gff line
        clean_line = line.strip().split('\t')

        # skip malformed lines (gff requires 9 columns)
        if len(clean_line) != 9:
            continue

        # extract feature type and count cds features
        feat_type = clean_line[2]
        if feat_type == 'CDS':
            cds_count += 1

        # only process cds features or transitions from cds
        if feat_type != 'CDS' and prev_feat_type != 'CDS':
            continue
        prev_feat_type = feat_type

        # initialize with first cds line
        if prev_line is None:
            prev_line = clean_line
            continue

        # extract gene information from previous line
        seqID, annotation_ID, parent, direction, gene_start, gene_end, nuc_length, aa_length = extract_gene_info(
            prev_line)

        # get next gene's parent and direction for cugo boundary analysis
        next_parent = clean_line[0]
        next_direction = clean_line[6]

        # determine cugo boundaries and classifications
        cugo, cugo_start, cugo_end, cugo_count, cugo_size_count = cugo_boundaries(
            direction, prev_direction, next_direction,
            parent, prev_parent, next_parent,
            cugo_count, cugo_size, cugo_size_count, prev_cugo
        )

        # store processed gene data (all 11 fields for database)
        reformat_data.append([seqID, parent, aa_length,
                              direction, annotation_ID, cugo])

        # update tracking variables for next iteration
        prev_line = clean_line
        prev_direction = direction
        prev_parent = parent
        prev_cugo = cugo

    # process the last gene in the file (no next gene to compare)
    if len(clean_line) == 9 and clean_line[2] == 'CDS':
        seqID, annotation_ID, parent, direction, gene_start, gene_end, nuc_length, aa_length = extract_gene_info(
            clean_line)

        # handle last gene - no next gene, so use none for next_parent
        cugo, cugo_start, cugo_end, cugo_count, cugo_size_count = cugo_boundaries(
            direction, prev_direction, None,  # no next direction for last gene
            parent, prev_parent, None,  # no next parent for last gene
            cugo_count, cugo_size, cugo_size_count, prev_cugo
        )

        # store final gene data (all 6 fields for database)
        reformat_data.append([seqID, parent, aa_length,
                              direction, annotation_ID, cugo])

        # finalize cugo size tracking
        cugo_size[cugo] = cugo_size_count

    # Return complete data for database storage
    return reformat_data

def parse_gff_for_update(gff_content: str) -> list:
    """
    Parses GFF content for UPDATE operations (seqID and annotation_ID only).

    Args:
        gff_content (str): GFF format file content as string

    Returns:
        list: List of tuples (annotation_ID, seqID)
    """
    update_data = []
    lines = gff_content.strip().split('\n')

    for line in lines:
        clean_line = line.strip().split('\t')

        if len(clean_line) != 9 or clean_line[2] != 'CDS':
            continue

        attributes = clean_line[8]

        # parse attributes
        attr_dict = {}
        for attr in attributes.split(';'):
            if '=' in attr:
                key, value = attr.split('=', 1)
                attr_dict[key] = value

        seqID = attr_dict.get('ID', '')
        annotation_ID = attr_dict.get('Name', '')

        if seqID and annotation_ID:
            update_data.append((annotation_ID, seqID))

    return update_data



def process_gff_file(filepath: str,
                     update_mode: bool = False):
    """
    Process a single GFF file (gzipped or plain) from disk and parse it.

    Args:
        filepath (str): Full path to the GFF file.

    Returns:
        list: Parsed gene records.
    """
    try:
        filepath = Path(filepath)

        # Open gzipped or plain text
        if filepath.suffix == ".gz":
            with gzip.open(filepath, "rt") as f:
                content = f.read()
        else:
            with open(filepath, "r", encoding="utf-8") as f:
                content = f.read()

        if update_mode:
            return parse_gff_for_update(content)
        else:
            return parse_single_gff(content)
    except Exception as e:
        # Optional: log the error
        # logger.warning(f"Failed to process {filepath}: {e}")
        return []

def process_tmhmm_file(tmhmm_filepath):
    try:
        with open(tmhmm_filepath, 'r') as f:
            tmhmm_data = []
            for line in f:
                if line.startswith('prot_ID'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    prot_id = parts[0]
                    try:
                        no_tmh = int(parts[1])
                    except ValueError:
                        no_tmh = 0
                    tmhmm_data.append((no_tmh, prot_id))
        return tmhmm_data
    except Exception:
        return []

def parse(cog_gff_tar_path: str,
          kegg_gff_tar_path: str = None,
          pfam_gff_tar_path: str = None,
          tmhmm_tar_path: str = None,
          output_dir: str = None,
          globdb_version: int = None,
          force: bool = False,
          ):
    """
    Main parsing function - processes COG first with full data, then updates with KEGG and Pfam.
    All annotations use Name= attribute from GFF files.
    """
    from logging import getLogger
    logger = getLogger(__name__)

    logger.info("Processing files: COG → TMHMM → KEGG → Pfam")

    # db_path = ensure_path(target=f"globdb_{globdb_version}_cugo.db")
    db_path = f"globdb_{globdb_version}_cugo.db"

    if output_dir is None:
        output_dir = '.'

    output_path = Path(output_dir)
    tempdir = output_path / 'tempdir'
    if tempdir.exists():
        shutil.rmtree(tempdir)
    tempdir.mkdir(parents=True, exist_ok=True)

    # Setup database
    conn = setup_database(db_path)

    # ===== STEP 1: Process COG GFF files =====
    logger.info("STEP 1: Processing COG GFF files...")
    subprocess.run(["tar", "-xzf", cog_gff_tar_path, "-C", str(tempdir)], check=True)

    cog_gff_files = list(tempdir.rglob("*.gff.gz")) + list(tempdir.rglob("*.gff"))
    logger.info(f"Found {len(cog_gff_files)} COG GFF files")

    for gff_file in tqdm(cog_gff_files, desc="Processing COG GFF"):
        gff_data = process_gff_file(str(gff_file), update_mode=False)
        if gff_data:
            conn.executemany("""
                INSERT OR REPLACE INTO all_data
                (seqID, parent_ID, aa_length, strand, COG_ID, cugo_number)
                VALUES (?, ?, ?, ?, ?, ?)
            """, gff_data)
            conn.commit()

    shutil.rmtree(tempdir)
    tempdir.mkdir(parents=True, exist_ok=True)

    # ===== STEP 2: Process TMHMM files =====
    if tmhmm_tar_path:
        logger.info("STEP 2: Processing TMHMM files...")
        subprocess.run(["tar", "-xzf", tmhmm_tar_path, "-C", str(tempdir)], check=True)

        tmhmm_files = list(tempdir.rglob("*_tmhmm_clean"))
        logger.info(f"Found {len(tmhmm_files)} TMHMM files")

        for tmhmm_file in tqdm(tmhmm_files, desc="Processing TMHMM"):
            tmhmm_data = process_tmhmm_file(str(tmhmm_file))
            if tmhmm_data:
                conn.executemany("""
                                 UPDATE all_data
                                 SET no_tmh = ?
                                 WHERE seqID = ?
                                 """, tmhmm_data)
                conn.commit()

        shutil.rmtree(tempdir)
        tempdir.mkdir(parents=True, exist_ok=True)

    # ===== STEP 3: Process KEGG GFF files =====
    if kegg_gff_tar_path:
        logger.info("STEP 3: Updating KEGG annotations...")
        subprocess.run(["tar", "-xzf", kegg_gff_tar_path, "-C", str(tempdir)], check=True)

        kegg_gff_files = list(tempdir.rglob("*.gff.gz")) + list(tempdir.rglob("*.gff"))
        logger.info(f"Found {len(kegg_gff_files)} KEGG GFF files")

        for gff_file in tqdm(kegg_gff_files, desc="Updating KEGG IDs"):
            kegg_data = process_gff_file(str(gff_file), update_mode=True)
            if kegg_data:
                conn.executemany("""
                                 UPDATE all_data
                                 SET KEGG_ID = ?
                                 WHERE seqID = ?
                                 """, kegg_data)
                conn.commit()

        shutil.rmtree(tempdir)
        tempdir.mkdir(parents=True, exist_ok=True)

    # ===== STEP 4: Process Pfam GFF files =====
    if pfam_gff_tar_path:
        logger.info("STEP 4: Updating Pfam annotations...")
        subprocess.run(["tar", "-xzf", pfam_gff_tar_path, "-C", str(tempdir)], check=True)

        pfam_gff_files = list(tempdir.rglob("*.gff.gz")) + list(tempdir.rglob("*.gff"))
        logger.info(f"Found {len(pfam_gff_files)} Pfam GFF files")

        for gff_file in tqdm(pfam_gff_files, desc="Updating Pfam IDs"):
            pfam_data = process_gff_file(str(gff_file), update_mode=True)
            if pfam_data:
                conn.executemany("""
                                 UPDATE all_data
                                 SET Pfam_ID = ?
                                 WHERE seqID = ?
                                 """, pfam_data)
                conn.commit()

        shutil.rmtree(tempdir)

    conn.close()
    logger.info("Database creation complete!")


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
    """
    Extracts top N most frequent values and their counts from each column.

    Args:
        df (pd.DataFrame): input dataframe to analyze
        top_n (int): number of top values to extract per column

    Returns:
        tuple: (id_df, count_df) where:
            - id_df (pd.DataFrame): top N values for each column
            - count_df (pd.DataFrame): corresponding counts for each value
    """
    top_ids = [[] for _ in range(top_n)]
    top_counts = [[] for _ in range(top_n)]

    for col in df.columns:
        # get value counts for current column
        counts = df[col].value_counts()
        for i in range(top_n):
            if i < len(counts):
                top_ids[i].append(counts.index[i])
                top_counts[i].append(counts.iloc[i])
            else:
                # fill with nan if fewer than top_n values exist
                top_ids[i].append(np.nan)
                top_counts[i].append(np.nan)

    id_df = pd.DataFrame(top_ids, columns=df.columns)
    count_df = pd.DataFrame(top_counts, columns=df.columns)
    return id_df, count_df


def plot_top_cogs_per_position(
        context_path: str,
        flank_lower: int,
        flank_upper: int,
        top_n: int,
        title: str = 'Top COGs per Position',
        save: bool = False,
        ax: plt.Axes = None,
        plot_path: str = None
):
    """
    Creates scatter plot showing top COG categories at each genomic position.

    Args:
        context_path (str): path to context data file
        flank_lower (int): lower bound of flanking region
        flank_upper (int): upper bound of flanking region
        top_n (int): number of top COGs to display per position
        title (str): plot title
        save (bool): whether to save plot to file
        ax (plt.Axes, optional): matplotlib axes object to plot on
        plot_path (str, optional): path to save plot file

    Returns:
        tuple: (positions, boundaries) where:
            - positions (list): list of position values
            - boundaries (list): list of position boundary values
    """
    # load context data and extract cog information
    cont = pd.read_csv(context_path, sep='\t', na_values='', keep_default_na=False, dtype=str)
    extract = cont.loc[cont['feat_type'] == 'COG_ID']

    # load cog color mapping
    script_dir = Path(__file__).resolve().parent
    color_yaml_path = script_dir / 'cog_colors.yaml'

    with open(color_yaml_path, 'r') as f:
        cog_color_map = yaml.safe_load(f)

    # select columns for specified flanking region
    flank_cols = [str(i) for i in range(flank_lower, flank_upper + 1)]
    cugo_df = extract[flank_cols]
    positions = [int(col) for col in cugo_df.columns]

    # extract top n values and counts for each position
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

    # calculate positioning for multiple ranks at each position
    subtick_width = 0.8 / top_n
    subtick_offset = (top_n - 1) * subtick_width / 2

    x_pos, y_values, cog_labels, point_colors = [], [], [], []
    all_xticks, all_xlabels = [], []

    # create figure if no axes provided
    if ax is None:
        figsize = (max(8, len(positions) * 0.8), 7)
        fig, ax = plt.subplots(figsize=figsize)

    # prepare data points for each position and rank
    for pos_idx, pos in enumerate(positions):
        for rank in range(top_n):
            cog_id = id_df.loc[rank, str(pos)]
            count = count_df.loc[rank, str(pos)]

            # calculate x position with offset for multiple ranks
            x = pos + (rank * subtick_width) - subtick_offset
            x_pos.append(x)
            y_values.append(count)
            all_xticks.append(x)

            if pd.isna(cog_id):
                cog_labels.append('nan')
                point_colors.append('#ffffff')
                all_xlabels.append('nan')
            else:
                cog_labels.append(cog_id)
                point_colors.append(cog_color_map.get(cog_id, '#aaaaaa'))
                all_xlabels.append(cog_id)

    # create scatter plot
    ax.scatter(
        x_pos,
        y_values,
        color=point_colors,
        s=150,
        zorder=3,
        edgecolors='black',
        linewidths=0.7
    )

    # set x-axis labels and ticks
    ax.set_xticks(all_xticks)
    ax.set_xticklabels(all_xlabels, rotation=90, fontsize=14, ha='center')

    # add vertical lines between positions
    for pos in positions[:-1]:
        boundary = pos + 0.5
        ax.axvline(boundary, color='gray', linestyle='--', linewidth=0.7, zorder=1)

    # add position labels at top
    max_y = max([y for y in y_values if not pd.isna(y)], default=0)
    label_offset = max_y * 0.05

    for pos in positions:
        ax.text(pos, max_y + label_offset, str(pos), ha='center', va='bottom',
                fontsize=16, fontweight='bold')

    # set axis limits and labels
    ax.set_ylim(-max_y * 0.05, max_y * 1.15)
    ax.set_xlim(positions[0] - 0.5, positions[-1] + 0.5)

    ax.tick_params(axis='y', labelsize=14)
    ax.set_ylabel('Count', fontsize=16)
    ax.set_title(title, fontsize=18, fontweight='bold')

    # save plot if requested
    if save:
        plt.tight_layout()
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')

    return positions, [pos + 0.5 for pos in positions[:-1]]


def plot_size_per_position(context_path: str,
                           flank_lower: int,
                           flank_upper: int,
                           title: str = 'aa_length Density per Position',
                           save: bool = False,
                           ax: Optional[plt.Axes] = None,
                           pos_boundaries: Optional[list] = None,
                           bin_width: int = 10,
                           y_range: int = None,
                           plot_path: str = None
                           ):
    """
    Creates 1D-density plot showing sequence length distribution across genomic positions.

    Args:
        context_path (str): path to context data file
        flank_lower (int): lower bound of flanking region
        flank_upper (int): upper bound of flanking region
        title (str): plot title
        save (bool): whether to save plot to file
        ax (plt.Axes, optional): matplotlib axes object to plot on
        pos_boundaries (list, optional): position boundaries to draw as vertical lines
        bin_width (int): width of length bins for histogram
        y_range (int, optional): maximum y-axis range to display
        plot_path (str, optional): path to save plot file

    Returns:
        None
    """
    # load context data and extract amino acid length information
    cont = pd.read_csv(context_path, sep='\t', na_values='', keep_default_na=False, dtype=str)
    extract = cont.loc[cont['feat_type'] == 'aa_length']

    # select columns for specified flanking region
    flank_cols = [str(i) for i in range(flank_lower, flank_upper + 1)]
    cugo_df = extract[flank_cols]
    positions = [int(col) for col in cugo_df.columns]

    # determine bin edges for length histograms
    all_lengths = cugo_df.values.flatten()
    all_lengths = all_lengths[~pd.isna(all_lengths)].astype(float)
    max_len = all_lengths.max()
    bin_edges = np.arange(0, max_len + bin_width, bin_width)
    n_bins = len(bin_edges) - 1

    # create histogram data for each position
    heat_data = np.zeros((n_bins, len(flank_cols)))
    position_counts = []

    for col_idx, col in enumerate(flank_cols):
        values = cugo_df[col].dropna().astype(float)
        hist, _ = np.histogram(values, bins=bin_edges)
        heat_data[:, col_idx] = hist
        position_counts.append(len(values))

    # create figure if no axes provided
    if ax is None:
        fig, ax = plt.subplots()

    # set up colormap and normalization
    cmap = get_cmap('Blues')
    norm = Normalize(vmin=0, vmax=heat_data.max())

    # draw heatmap rectangles
    rect_width = 0.8
    for col_idx, pos in enumerate(positions):
        for y in range(heat_data.shape[0]):
            val = heat_data[y, col_idx]
            color = cmap(norm(val))
            ax.add_patch(plt.Rectangle(
                (pos - rect_width / 2, y), rect_width, 1, color=color, linewidth=0
            ))

    # add position boundaries if provided
    if pos_boundaries is not None:
        for boundary in pos_boundaries:
            ax.axvline(boundary, color='gray', linestyle='--', linewidth=0.7, zorder=1)

    # set x-axis labels with position and count information
    ax.set_xlim(positions[0] - 0.5, positions[-1] + 0.5)
    xtick_labels = [f'{pos}\n(n={count})' for pos, count in zip(positions, position_counts)]
    ax.set_xticks(positions)
    ax.set_xticklabels(xtick_labels, fontsize=14)

    # set y-axis range
    if y_range:
        n_bins = int(y_range / bin_width)

    ax.set_ylim(0, n_bins)

    # Set y-ticks at 0, middle, and maximum
    max_length = n_bins * bin_width
    ax.set_yticks([0, n_bins/2, n_bins])
    ax.set_yticklabels(['0', f'{int(max_length/2)}', f'{int(max_length)}'], fontsize=16)

    # set axis labels and title
    ax.set_xlabel('Position', fontsize=18)
    ax.set_ylabel(f'length (bin size: {bin_width})', fontsize=18)
    ax.set_title(title, fontsize=20)

    # create colorbar
    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])

    plt.tight_layout()
    if save:
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')

def plot_tmh_per_position(context_path: str,
                          flank_lower: int,
                          flank_upper: int,
                          title: str = 'no_TMH Distribution per Position',
                          save: bool = False,
                          ax: Optional[plt.Axes] = None,
                          pos_boundaries: Optional[list] = None,
                          y_range: int = None,
                          plot_path: str = None
                          ):
    """
    Creates 1D-density plot showing transmembrane helix count distribution across genomic positions.

    Args:
        context_path (str): path to context data file
        flank_lower (int): lower bound of flanking region
        flank_upper (int): upper bound of flanking region
        title (str): plot title
        save (bool): whether to save plot to file
        ax (plt.Axes, optional): matplotlib axes object to plot on
        pos_boundaries (list, optional): position boundaries to draw as vertical lines
        y_range (int, optional): maximum y-axis range to display
        plot_path (str, optional): path to save plot file

    Returns:
        None
    """
    # load context data and extract no_TMH information
    cont = pd.read_csv(context_path, sep='\t', na_values='', keep_default_na=False, dtype=str)
    extract = cont.loc[cont['feat_type'] == 'no_TMH']

    # select columns for specified flanking region
    flank_cols = [str(i) for i in range(flank_lower, flank_upper + 1)]
    cugo_df = extract[flank_cols]
    positions = [int(col) for col in cugo_df.columns]

    # determine bin edges for TMH counts (integers)
    all_tmh = cugo_df.values.flatten()
    all_tmh = all_tmh[~pd.isna(all_tmh)].astype(float)
    max_tmh = int(all_tmh.max())
    bin_edges = np.arange(0, max_tmh + 2, 1)  # bins for 0, 1, 2, ..., max_tmh
    n_bins = len(bin_edges) - 1

    # create histogram data for each position
    heat_data = np.zeros((n_bins, len(flank_cols)))
    position_counts = []

    for col_idx, col in enumerate(flank_cols):
        values = cugo_df[col].dropna().astype(float)
        hist, _ = np.histogram(values, bins=bin_edges)
        heat_data[:, col_idx] = hist
        position_counts.append(len(values))

    # create figure if no axes provided
    if ax is None:
        fig, ax = plt.subplots()

    # set up colormap and normalization
    cmap = get_cmap('Reds')
    norm = Normalize(vmin=0, vmax=heat_data.max())

    # draw heatmap rectangles
    rect_width = 0.8
    for col_idx, pos in enumerate(positions):
        for y in range(heat_data.shape[0]):
            val = heat_data[y, col_idx]
            color = cmap(norm(val))
            ax.add_patch(plt.Rectangle(
                (pos - rect_width / 2, y), rect_width, 1, color=color, linewidth=0
            ))

    # add position boundaries if provided
    if pos_boundaries is not None:
        for boundary in pos_boundaries:
            ax.axvline(boundary, color='gray', linestyle='--', linewidth=0.7, zorder=1)

    # set x-axis labels with position and count information
    ax.set_xlim(positions[0] - 0.5, positions[-1] + 0.5)
    xtick_labels = [f'{pos}\n(n={count})' for pos, count in zip(positions, position_counts)]
    ax.set_xticks(positions)
    ax.set_xticklabels(xtick_labels, fontsize=14)

    # set y-axis range
    if y_range:
        n_bins = min(y_range, n_bins)

    ax.set_ylim(0, n_bins)

    # Set y-ticks at 0, middle, and maximum (showing actual TMH values)
    ax.set_yticks([0, n_bins/2, n_bins])
    ax.set_yticklabels(['0', f'{int(n_bins/2)}', f'{int(n_bins)}'], fontsize=16)

    # set axis labels and title
    ax.set_xlabel('Position', fontsize=18)
    ax.set_ylabel('no_TMH', fontsize=18)
    ax.set_title(title, fontsize=20)

    # create colorbar
    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])

    plt.tight_layout()
    if save:
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')

# ======================================
# CUGO COMMAND LINE TOOLS
# ======================================
def process_target_context(target_id, parent_id, target_cugo, strand, cugo_range,
                           context_data, headers, seq_idx):
    """Build context window for a single target protein."""
    if parent_id not in context_data:
        return []

    # Get rows in CUGO number order
    parent_rows = sorted(context_data[parent_id], key=lambda x: int(x[headers.index('CUGO_number')]))

    # Filter to exact range
    context_window = []
    for row in parent_rows:
        cugo_num = int(row[headers.index('CUGO_number')])
        if target_cugo - cugo_range <= cugo_num <= target_cugo + cugo_range:
            context_window.append(row)

    if not context_window:
        return []

    # Reverse for negative strand
    if strand == '-':
        context_window = context_window[::-1]

    # Find target index
    target_index = None
    for i, row in enumerate(context_window):
        if row[seq_idx] == target_id:
            target_index = i
            break

    if target_index is None:
        return []

    # Transpose: each column becomes a row
    results = []
    for col_idx, header in enumerate(headers):
        row_data = {'feat_type': header}

        for i, row in enumerate(context_window):
            centered_idx = i - target_index
            row_data[str(centered_idx)] = row[col_idx]

        results.append(row_data)

    return results


def write_context_output(results, output_file):
    """Write the transposed results to output file."""
    if not results:
        return

    # Get all column indices
    all_indices = set()
    for row in results:
        for key in row.keys():
            if key != 'feat_type':
                all_indices.add(int(key))

    sorted_indices = sorted(all_indices)

    with open(output_file, 'w') as f:
        # Header
        f.write('feat_type\t' + '\t'.join(str(idx) for idx in sorted_indices) + '\n')

        # Data rows
        for row in results:
            line_parts = [row['feat_type']]
            for idx in sorted_indices:
                line_parts.append(row.get(str(idx), ''))
            f.write('\t'.join(map(str, line_parts)) + '\n')


def fetch_seqid_batch(batch, cugo_path):
    conn = sqlite3.connect(cugo_path)
    cursor = conn.cursor()
    placeholders = ",".join("?" * len(batch))
    query = f"""
        SELECT seqID, parent_ID, aa_length, strand, COG_ID, CUGO_number, no_TMH
        FROM all_data
        WHERE seqID IN ({placeholders})
    """
    cursor.execute(query, batch)
    rows = cursor.fetchall()
    conn.close()
    return rows

def fetch_parent_context(parent_ID, needed_numbers, cugo_path, batch_size=500):
    conn = sqlite3.connect(cugo_path)
    cursor = conn.cursor()
    rows_out = []
    needed_list = list(needed_numbers)
    for j in range(0, len(needed_list), batch_size):
        batch = needed_list[j:j + batch_size]
        placeholders = ",".join("?" * len(batch))
        query = f"""
            SELECT seqID, parent_ID, aa_length, strand, COG_ID, CUGO_number, no_TMH
            FROM all_data
            WHERE parent_ID = ? AND CUGO_number IN ({placeholders})
            ORDER BY CUGO_number
        """
        cursor.execute(query, (parent_ID, *batch))
        rows_out.extend(cursor.fetchall())
    conn.close()
    return parent_ID, rows_out

def context(fasta_path: str,
            cugo_path: str,
            cugo_range: int,
            output_dir: str,
            threads: int = 1,
            force: bool = False,
            ):
    if fasta_path:
        protein_name = determine_dataset_name(fasta_path, '.', 0)
        output_file = ensure_path(output_dir, f'{protein_name}_context.tsv', force=force)
        sequences = read_fasta_to_dict(fasta_path)
        protein_identifiers = set(sequences.keys())  # O(1) lookups

        log_file = ensure_path(output_dir, f'{protein_name}_missing_files.log', force=force)
        logging.basicConfig(
            filename=log_file,
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
    else:
        raise ValueError('You must provide a FASTA file.')

    BATCH_SIZE = 500
    target_rows = []
    protein_list = list(protein_identifiers)
    batches = [protein_list[i:i + BATCH_SIZE] for i in range(0, len(protein_list), BATCH_SIZE)]

    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(fetch_seqid_batch, batch, cugo_path): batch for batch in batches}
        for future in tqdm(as_completed(futures), total=len(futures), desc="Fetching seqIDs"):
            target_rows.extend(future.result())

    if not target_rows:
        logging.warning("No target proteins found in the database.")
        return None

    target_contexts = {}
    context_ranges = defaultdict(set)

    for seqID, parent_ID, aa_length, strand, COG_ID, CUGO_number, no_TMH in target_rows:
        target_contexts[seqID] = (parent_ID, CUGO_number, strand)
        for i in range(CUGO_number - cugo_range, CUGO_number + cugo_range + 1):
            context_ranges[parent_ID].add(i)

    headers = ["seqID", "parent_ID", "aa_length", "strand", "COG_ID", "CUGO_number", "no_TMH"]
    context_data = defaultdict(list)

    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {
            executor.submit(fetch_parent_context, pid, numbers, cugo_path, BATCH_SIZE): pid
            for pid, numbers in context_ranges.items()
        }
        for future in tqdm(as_completed(futures), total=len(futures), desc="Fetching parent contexts"):
            parent_ID, rows = future.result()
            if rows:
                context_data[parent_ID].extend(rows)

    all_results = []
    for seqID, (parent_ID, target_cugo, strand) in target_contexts.items():
        results = process_target_context(
            target_id=seqID,
            parent_id=parent_ID,
            target_cugo=target_cugo,
            strand=strand,
            cugo_range=cugo_range,
            context_data=context_data,
            headers=headers,
            seq_idx=0,
        )
        all_results.extend(results)

    if all_results:
        write_context_output(all_results, output_file)
        logging.info(f"Context written to {output_file}")
        return output_file
    else:
        logging.info("No context data found.")
        return None


def cugo_plot(context_path: str,
              flank_lower: int,
              flank_upper: int,
              top_n: int,
              output: str,
              cugo: bool = False,
              size: bool = False,
              tmh: bool = False,
              all_plots: bool = False,
              bin_width: int = 10,
              y_range: int = None,
              tmh_y_range: int = None,
              force: bool = False
              ):
    """
    Generate plots for genomic context analysis: COG distribution, protein size, and TMH counts.

    Args:
        context_path: Path to context TSV file
        flank_lower: Lower flank boundary for plotting
        flank_upper: Upper flank boundary for plotting
        top_n: Number of top COGs to display
        output: Output directory for plots
        cugo: Whether to generate COG-only plot
        size: Whether to generate size-only plot
        tmh: Whether to generate TMH-only plot
        all_plots: Whether to generate combined plot
        bin_width: Bin width for size plots
        y_range: Y-axis range for size plots
        tmh_y_range: Y-axis range for TMH plots
        force: Whether to overwrite existing files
    """
    dataset_name = context_path.removesuffix('_context.tsv')

    # generate cog-only plot
    if cugo:
        if not top_n:
            logger.error('top_n is required if cugo is True')
        else:
            cugo_plot_path = ensure_path(output, f'{dataset_name}_cugo_only.png', force=force)
            plot_top_cogs_per_position(context_path=context_path, flank_lower=flank_lower,
                                       flank_upper=flank_upper, top_n=top_n, save=True,
                                       plot_path=cugo_plot_path)
            logger.INFO(f"Plot saved to {cugo_plot_path}")


# generate size-only plot
    if size:
        size_plot_path = ensure_path(output, f'{dataset_name}_size_only.png', force=force)
        plot_size_per_position(context_path=context_path, flank_lower=flank_lower,
                               flank_upper=flank_upper, save=True, bin_width=bin_width,
                               plot_path=size_plot_path, y_range=y_range)
        logger.INFO(f"Plot saved to {size_plot_path}")


# generate tmh-only plot
    if tmh:
        tmh_plot_path = ensure_path(output, f'{dataset_name}_tmh_only.png', force=force)
        plot_tmh_per_position(context_path=context_path, flank_lower=flank_lower,
                              flank_upper=flank_upper, save=True, y_range=tmh_y_range,
                              plot_path=tmh_plot_path)
        logger.INFO(f"Plot saved to {tmh_plot_path}")


# generate combined plot
    if all_plots:
        if not top_n:
            logger.error('top_n is required if all_plots is True')
        else:
            all_plot_path = ensure_path(output, f'{dataset_name}_cugo.png', force=force)

            # calculate dynamic figure width
            width = max(8, int((flank_upper - flank_lower + 1) * top_n * 0.6))
            figsize = (width, 16)
            fig, (ax1, ax2, ax3) = plt.subplots(3, 1,
                                                figsize=figsize,
                                                gridspec_kw={'height_ratios': [3, 2, 1.5]})

            # plot top cogs per position
            pos_centers, pos_boundaries = plot_top_cogs_per_position(
                context_path=context_path,
                flank_lower=flank_lower,
                flank_upper=flank_upper,
                top_n=top_n,
                title='Top COGs per Position',
                ax=ax1,
            )

            # plot protein size distribution
            plot_size_per_position(
                context_path=context_path,
                flank_lower=flank_lower,
                flank_upper=flank_upper,
                title='aa_length Density per Position',
                ax=ax2,
                pos_boundaries=pos_boundaries,
                bin_width=bin_width,
                y_range=y_range
            )

            # plot TMH distribution
            plot_tmh_per_position(
                context_path=context_path,
                flank_lower=flank_lower,
                flank_upper=flank_upper,
                title='no_TMH Distribution per Position',
                ax=ax3,
                pos_boundaries=pos_boundaries,
                y_range=tmh_y_range
            )

            ax1.tick_params(labelbottom=True)

            plt.tight_layout()
            plt.subplots_adjust(hspace=0.4)
            plt.savefig(all_plot_path, dpi=300, bbox_inches='tight')
            logger.info(f"Plot saved to {all_plot_path}")

# ======================================
# CUGO COMMAND LINE WORKFLOW
# ======================================

def cugo(cugo_path: str,
         cugo_range: int,
         fasta_path: str,
         output_dir: str,
         flank_lower: int,
         flank_upper: int,
         top_n: int,
         threads: int = 1,
         force: bool = False,
         bin_width: int = 10,
         y_range: int = None,
         tmh_y_range: int = None):
    """
    Complete CUGO workflow: generate context data and create comprehensive plots.

    Args:
        cugo_path: Path to CUGO file
        cugo_range: Range around target protein for context
        output_dir: Directory for output files
        flank_lower: Lower flank boundary for plotting
        flank_upper: Upper flank boundary for plotting
        top_n: Number of top COGs to display
        force: Whether to overwrite existing files
        fasta_path: Optional path to FASTA file
        bin_width: Bin width for size plots
        y_range: Y-axis range for size plots

    Returns:
        tuple: (context_file_path, plot_file_path)
    """
    # generate context data
    context_file = context(
        fasta_path=fasta_path,
        cugo_path=cugo_path,
        cugo_range=cugo_range,
        output_dir=output_dir,
        threads=threads,
        force=force
    )

    if context_file is None:
        logger.error("Context generation failed - no plots will be created")
        return None, None

    # create comprehensive plots
    cugo_plot(
        context_path=context_file,
        flank_lower=flank_lower,
        flank_upper=flank_upper,
        top_n=top_n,
        output=output_dir,
        all_plots=True,
        bin_width=bin_width,
        y_range=y_range,
        tmh_y_range=tmh_y_range,
        force=force
    )

    # determine plot file path
    dataset_name = context_file.removesuffix('_context.tsv')
    plot_file = f"{dataset_name}_cugo.png"

    return context_file, plot_file
