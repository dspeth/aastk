#!/usr/bin/env python3

from .util import *
from .database import ANNOTATION_COLUMNS

import pandas as pd
import logging
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable, get_cmap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from typing import Optional
from pathlib import Path
import yaml
import gzip
import sqlite3
import hashlib
from tqdm import tqdm
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed


logger = logging.getLogger(__name__)

# ======================================
# CUGO context functions and CLI tool
# ======================================
def process_target_context(target_id, parent_id, target_cugo, strand, cugo_range,
                           context_data, headers, seq_idx):
    """Build context window for a single target protein."""
    if parent_id not in context_data:
        return []

    parent_rows = sorted(context_data[parent_id], key=lambda x: int(x[headers.index('CUGO_number')]))

    context_window = []
    for row in parent_rows:
        cugo_num = int(row[headers.index('CUGO_number')])
        if target_cugo - cugo_range <= cugo_num <= target_cugo + cugo_range:
            context_window.append(row)

    if not context_window:
        return []

    if strand == '-':
        context_window = context_window[::-1]

    target_index = None
    for i, row in enumerate(context_window):
        if row[seq_idx] == target_id:
            target_index = i
            break

    if target_index is None:
        return []

    # build rows with target_id and position
    results = []
    for i, row in enumerate(context_window):
        position = i - target_index
        row_dict = {'target_id': target_id, 'position': position}
        for col_idx, header in enumerate(headers):
            row_dict[header] = row[col_idx]
        results.append(row_dict)

    return results


def write_context_output(all_results, annotation, output_file):
    """Write all results to one TSV."""
    if not all_results:
        return

    headers = ['target_id', 'position', 'seqID', 'parent_ID', 'aa_length',
               'strand', f'{annotation}', 'CUGO_number', 'no_TMH']

    with open(output_file, 'w') as f:
        f.write('\t'.join(headers) + '\n')
        for row in all_results:
            line = [str(row.get(h, '')) for h in headers]
            f.write('\t'.join(line) + '\n')


def fetch_seqid_batch(batch, annotation, db_path):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    placeholders = ",".join("?" * len(batch))
    query = f"""
        SELECT seqID, parent_ID, aa_length, strand, {annotation}, CUGO_number, no_TMH
        FROM protein_data
        WHERE seqID IN ({placeholders})
    """
    cursor.execute(query, batch)
    rows = cursor.fetchall()
    conn.close()
    return rows

def fetch_parent_context(parent_ID, needed_numbers, annotation, db_path, batch_size=500):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    rows_out = []
    needed_list = list(needed_numbers)
    for j in range(0, len(needed_list), batch_size):
        batch = needed_list[j:j + batch_size]
        placeholders = ",".join("?" * len(batch))
        query = f"""
            SELECT seqID, parent_ID, aa_length, strand, {annotation}, CUGO_number, no_TMH
            FROM protein_data
            WHERE parent_ID = ? AND CUGO_number IN ({placeholders})
            ORDER BY CUGO_number
        """
        cursor.execute(query, (parent_ID, *batch))
        rows_out.extend(cursor.fetchall())
    conn.close()
    return parent_ID, rows_out

def context(fasta: str,
            id_list: str,
            db_path: str,
            cugo_range: int,
            output_dir: str,
            annotation: str = 'COG_ID',
            threads: int = 1,
            force: bool = False,
            ):
    if annotation not in ANNOTATION_COLUMNS:
        raise ValueError(f'Invalid annotation. Please select one of the following annotations: {",".join(annotation_columns)}')

    if fasta:
        protein_name = determine_dataset_name(fasta, '.', 0)
        output_file = ensure_path(output_dir, f'{protein_name}_{annotation}_context.tsv', force=force)
        sequences = read_fasta_to_dict(fasta)
        protein_identifiers = set(sequences.keys())
    elif id_list:
        protein_name = determine_dataset_name(id_list, '.', 0)
        output_file = ensure_path(output_dir, f'{protein_name}_{annotation}_context.tsv', force=force)
        with open(id_list, 'r') as f:
            protein_identifiers = [line.strip() for line in f]
    else:
        raise ValueError('You must provide a FASTA file.')

    BATCH_SIZE = 500
    target_rows = []
    protein_list = list(protein_identifiers)
    batches = [protein_list[i:i + BATCH_SIZE] for i in range(0, len(protein_list), BATCH_SIZE)]

    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(fetch_seqid_batch, batch, annotation, db_path): batch for batch in batches}
        for future in tqdm(as_completed(futures), total=len(futures), desc="Fetching seqIDs"):
            target_rows.extend(future.result())

    if not target_rows:
        logging.warning("No target proteins found in the database.")
        return None

    target_contexts = {}
    context_ranges = defaultdict(set)

    for seqID, parent_ID, aa_length, strand, annotation_id, CUGO_number, no_TMH in target_rows:
        target_contexts[seqID] = (parent_ID, CUGO_number, strand)
        for i in range(CUGO_number - cugo_range, CUGO_number + cugo_range + 1):
            context_ranges[parent_ID].add(i)

    headers = ["seqID", "parent_ID", "aa_length", "strand", f"{annotation}", "CUGO_number", "no_TMH"]
    context_data = defaultdict(list)

    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {
            executor.submit(fetch_parent_context, pid, numbers, annotation, db_path, BATCH_SIZE): pid
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
        write_context_output(all_results, annotation, output_file)
        logging.info(f"Context written to {output_file}")
        return output_file
    else:
        logging.info("No context data found.")
        return None



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


def plot_top_annotations_per_position(
        context_path: str,
        flank_lower: int,
        flank_upper: int,
        top_n: int,
        annotation: str,
        save: bool = False,
        ax: plt.Axes = None,
        plot_path: str = None
):
    """
    Creates scatter plot showing top annotation categories at each genomic position.
    """
    title = f'Top {annotation}s per position'

    # Load context data
    cont = pd.read_csv(context_path, sep='\t', dtype=str)
    cont['position'] = pd.to_numeric(cont['position'], errors='coerce')
    cont = cont[(cont['position'] >= flank_lower) & (cont['position'] <= flank_upper)]

    positions = sorted(cont['position'].unique())

    # Load annotation color mapping
    script_dir = Path(__file__).resolve().parent
    yaml_dir = script_dir / "yaml"
    annotation_lower = annotation.lower().replace('_id', '')
    color_yaml_path = yaml_dir / f'{annotation_lower}_colors.yaml'

    # dynamically color by first 6 digits of md5 hash
    unique_annotations = cont[f'{annotation}'].dropna().unique()
    annotation_color_map = {
        ann_id: '#' + hashlib.md5(str(ann_id).encode()).hexdigest()[:6]
        for i, ann_id in enumerate(unique_annotations)
    }

    # Extract top N annotations per position
    top_ids = [[] for _ in range(top_n)]
    top_counts = [[] for _ in range(top_n)]

    for pos in positions:
        pos_data = cont[cont['position'] == pos]
        counts = pos_data[f'{annotation}'].value_counts(dropna=False)  # include NA
        for i in range(top_n):
            if i < len(counts):
                annotation_id = counts.index[i]
                top_ids[i].append(annotation_id)
                top_counts[i].append(counts.iloc[i])
            else:
                # Only append nothing if no real data in this rank
                top_ids[i].append(None)
                top_counts[i].append(None)

    # X-position offsets for multiple ranks
    subtick_width = 0.8 / top_n
    subtick_offset = (top_n - 1) * subtick_width / 2

    x_pos, y_values, annotation_labels, point_colors = [], [], [], []
    all_xticks, all_xlabels = [], []

    if ax is None:
        figsize = (max(8, len(positions) * 0.8), 7)
        fig, ax = plt.subplots(figsize=figsize)

    for pos_idx, pos in enumerate(positions):
        for rank in range(top_n):
            annotation_id = top_ids[rank][pos_idx]
            count = top_counts[rank][pos_idx]

            if annotation_id is None or count is None:
                continue  # skip completely empty ranks

            x = pos + (rank * subtick_width) - subtick_offset
            x_pos.append(x)
            y_values.append(count)
            all_xticks.append(x)

            if annotation_id == 'NA':
                annotation_labels.append('NA')
                point_colors.append('#cccccc')  # gray for NA
                all_xlabels.append('NA')
            else:
                annotation_labels.append(annotation_id)
                point_colors.append(annotation_color_map.get(annotation_id, '#aaaaaa'))
                all_xlabels.append(annotation_id)

    # Scatter plot
    ax.scatter(
        x_pos,
        y_values,
        color=point_colors,
        s=300,
        zorder=3,
        edgecolors='black',
        linewidths=0.7
    )

    # X-axis labels
    ax.set_xticks(all_xticks)
    ax.set_xticklabels(all_xlabels, rotation=90, fontsize=14, ha='center')

    # Vertical lines between positions
    for pos in positions[:-1]:
        ax.axvline(pos + 0.5, color='gray', linestyle='--', linewidth=0.7, zorder=1)

    # Position labels at top
    max_y = max([y for y in y_values if y is not None], default=0)
    label_offset = max_y * 0.05
    for pos in positions:
        ax.text(pos, max_y + label_offset, str(pos), ha='center', va='bottom',
                fontsize=16, fontweight='bold')

    # Axis limits and labels
    ax.set_ylim(-max_y * 0.05, max_y * 1.15)
    ax.set_xlim(positions[0] - 0.5, positions[-1] + 0.5)
    ax.set_xlabel('position', fontsize=18)
    ax.set_ylabel('Count', fontsize=18)
    ax.tick_params(axis='y', labelsize=16)
    ax.set_title(title, fontsize=20)

    # Save if requested
    if save and plot_path is not None:
        #plt.tight_layout()
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')

    return positions, [pos + 0.5 for pos in positions[:-1]]


def plot_size_per_position(context_path: str,
                           flank_lower: int,
                           flank_upper: int,
                           title: str = 'Density of amino acid length per position',
                           save: bool = False,
                           ax: Optional[plt.Axes] = None,
                           pos_boundaries: Optional[list] = None,
                           bin_width: int = 50,
                           y_range: int = None,
                           plot_path: str = None
                           ):
    """
    Creates 1D-density plot showing sequence length distribution across genomic positions.
    """
    # Load context data in long format
    cont = pd.read_csv(context_path, sep='\t')

    # Filter to position range
    cont = cont[(cont['position'] >= flank_lower) & (cont['position'] <= flank_upper)]

    # Convert aa_length to numeric
    cont['aa_length'] = pd.to_numeric(cont['aa_length'], errors='coerce')

    positions = sorted(cont['position'].unique())

    # Determine bin edges for length histograms
    all_lengths = cont['aa_length'].dropna()
    max_len = all_lengths.max()
    bin_edges = np.arange(0, max_len + bin_width, bin_width)
    n_bins = len(bin_edges) - 1

    # Create histogram data for each position
    heat_data = np.zeros((n_bins, len(positions)))
    position_counts = []

    for col_idx, pos in enumerate(positions):
        values = cont[cont['position'] == pos]['aa_length'].dropna()
        hist, _ = np.histogram(values, bins=bin_edges)
        heat_data[:, col_idx] = hist
        position_counts.append(len(values))

    # Create figure if no axes provided
    if ax is None:
        fig, ax = plt.subplots()

    # Set up colormap and normalization
    cmap = get_cmap('Blues')
    norm = Normalize(vmin=0, vmax=heat_data.max())

    # Draw heatmap rectangles
    rect_width = 0.8
    for col_idx, pos in enumerate(positions):
        for y in range(heat_data.shape[0]):
            val = heat_data[y, col_idx]
            color = cmap(norm(val))
            ax.add_patch(plt.Rectangle(
                (pos - rect_width / 2, y), rect_width, 1, color=color, linewidth=0
            ))

    # Add position boundaries if provided
    if pos_boundaries is not None:
        for boundary in pos_boundaries:
            ax.axvline(boundary, color='gray', linestyle='--', linewidth=0.7, zorder=1)

    # Set x-axis labels with position and count information
    ax.set_xlim(positions[0] - 0.5, positions[-1] + 0.5)
    xtick_labels = [f'{pos}\n(n={count})' for pos, count in zip(positions, position_counts)]
    ax.set_xticks(positions)
    ax.set_xticklabels(xtick_labels, fontsize=14)

    # Set y-axis range
    if y_range:
        n_bins = int(y_range / bin_width)

    ax.set_ylim(0, n_bins)

    # Set y-ticks at 0, middle, and maximum
    max_length = n_bins * bin_width
    ax.set_yticks([0, n_bins/2, n_bins])
    ax.set_yticklabels(['0', f'{int(max_length/2)}', f'{int(max_length)}'], fontsize=16)

    # Set axis labels and title
    ax.set_xlabel('position', fontsize=18)
    ax.set_ylabel(f'Length (bin size: {bin_width})', fontsize=18)
    ax.set_title(title, fontsize=20)

    if save:
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    return norm, cmap


def plot_tmh_per_position(context_path: str,
                          flank_lower: int,
                          flank_upper: int,
                          title: str = 'Density of transmembrane helix count per position',
                          save: bool = False,
                          ax: Optional[plt.Axes] = None,
                          pos_boundaries: Optional[list] = None,
                          y_range: int = None,
                          plot_path: str = None
                          ):
    """
    Creates 1D-density plot showing transmembrane helix count distribution across genomic positions.
    """
    # Load context data in long format
    cont = pd.read_csv(context_path, sep='\t')

    # Filter to position range
    cont = cont[(cont['position'] >= flank_lower) & (cont['position'] <= flank_upper)]

    # Convert no_TMH to numeric
    cont['no_TMH'] = pd.to_numeric(cont['no_TMH'], errors='coerce')

    positions = sorted(cont['position'].unique())

    # Determine bin edges for TMH counts (integers)
    all_tmh = cont['no_TMH'].dropna()
    max_tmh = int(all_tmh.max())
    bin_edges = np.arange(0, max_tmh + 2, 1)  # bins for 0, 1, 2, ..., max_tmh
    n_bins = len(bin_edges) - 1

    # Create histogram data for each position
    heat_data = np.zeros((n_bins, len(positions)))
    position_counts = []

    for col_idx, pos in enumerate(positions):
        values = cont[cont['position'] == pos]['no_TMH'].dropna()
        hist, _ = np.histogram(values, bins=bin_edges)
        heat_data[:, col_idx] = hist
        position_counts.append(len(values))

    # Create figure if no axes provided
    if ax is None:
        fig, ax = plt.subplots()

    # Set up colormap and normalization
    cmap = get_cmap('Reds')
    norm = Normalize(vmin=0, vmax=heat_data.max())

    # Draw heatmap rectangles
    rect_width = 0.8
    for col_idx, pos in enumerate(positions):
        for y in range(heat_data.shape[0]):
            val = heat_data[y, col_idx]
            color = cmap(norm(val))
            ax.add_patch(plt.Rectangle(
                (pos - rect_width / 2, y), rect_width, 1, color=color, linewidth=0
            ))

    # Add position boundaries if provided
    if pos_boundaries is not None:
        for boundary in pos_boundaries:
            ax.axvline(boundary, color='gray', linestyle='--', linewidth=0.7, zorder=1)

    # Set x-axis labels with position and count information
    ax.set_xlim(positions[0] - 0.5, positions[-1] + 0.5)
    xtick_labels = [f'{pos}\n(n={count})' for pos, count in zip(positions, position_counts)]
    ax.set_xticks(positions)
    ax.set_xticklabels(xtick_labels, fontsize=14)

    # Set y-axis range
    if y_range:
        n_bins = min(y_range, n_bins)

    ax.set_ylim(0, n_bins)

    # Set y-ticks at 0, middle, and maximum (showing actual TMH values)
    ax.set_yticks([0, n_bins/2, n_bins])
    ax.set_yticklabels(['0', f'{int(n_bins/2)}', f'{int(n_bins)}'], fontsize=16)

    # Set axis labels and title
    ax.set_xlabel('position', fontsize=18)
    ax.set_ylabel('Transmembrane helix count', fontsize=18)
    ax.set_title(title, fontsize=20)

    if save:
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')

    return norm, cmap

def cugo_plot(context_path: str,
              flank_lower: int,
              flank_upper: int,
              output: str,
              annotation: str = 'COG_ID',
              top_n: int = 3,
              cugo: bool = False,
              size: bool = False,
              tmh: bool = False,
              all_plots: bool = False,
              bin_width: int = 10,
              y_range: int = None,
              tmh_y_range: int = None,
              svg: bool = False,
              force: bool = False
              ):
    """
    Generate plots for genomic context analysis: annotation distribution, protein size, and TMH counts.

    Args:
        context_path: Path to context TSV file
        flank_lower: Lower flank boundary for plotting
        flank_upper: Upper flank boundary for plotting
        output: Output directory for plots
        annotation: Selects annotation type for analysis (default: COG_ID)
        top_n: Number of top annotation to display (default: 3)
        cugo: Whether to generate annotation-only plot
        size: Whether to generate size-only plot
        tmh: Whether to generate TMH-only plot
        all_plots: Whether to generate combined plot
        bin_width: Bin width for size plots
        y_range: Y-axis range for size plots
        tmh_y_range: Y-axis range for TMH plots
        svg: Generate plot in SVG format
        force: Whether to overwrite existing files
    """
    if annotation not in ANNOTATION_COLUMNS:
        raise ValueError(f'Invalid annotation. Please select one of the following annotations: {",".join(annotation_columns)}')

    dataset_name = determine_dataset_name(context_path, '.', 0, '_context')

    # generate annotation-only plot
    if cugo:
        logger.info(f'Plotting top {top_n} annotations per position.')
        if svg:
            cugo_plot_path = ensure_path(output, f'{dataset_name}_cugo_only.svg', force=force)
        else:
            cugo_plot_path = ensure_path(output, f'{dataset_name}_cugo_only.png', force=force)

        plot_top_annotations_per_position(context_path=context_path, flank_lower=flank_lower,
                                   flank_upper=flank_upper, annotation=annotation, top_n=top_n, save=True,
                                   plot_path=cugo_plot_path)
        (logger.info
         (f"Plot saved to {cugo_plot_path}"))


    # generate size-only plot
    if size:
        if svg:
            size_plot_path = ensure_path(output, f'{dataset_name}_size_only.svg', force=force)
        else:
            size_plot_path = ensure_path(output, f'{dataset_name}_size_only.png', force=force)

        norm_size, cmap_size = plot_size_per_position(context_path=context_path, flank_lower=flank_lower,
                                                      flank_upper=flank_upper, save=True, bin_width=bin_width,
                                                      plot_path=size_plot_path, y_range=y_range)
        logger.info(f"Plot saved to {size_plot_path}")


    # generate tmh-only plot
    if tmh:
        if svg:
            tmh_plot_path = ensure_path(output, f'{dataset_name}_tmh_only.svg', force=force)
        else:
            tmh_plot_path = ensure_path(output, f'{dataset_name}_tmh_only.png', force=force)

        norm_tmh, cmap_tmh = plot_tmh_per_position(context_path=context_path, flank_lower=flank_lower,
                                                   flank_upper=flank_upper, save=True, y_range=tmh_y_range,
                                                   plot_path=tmh_plot_path)
        logger.info(f"Plot saved to {tmh_plot_path}")


    # generate combined plot
    if all_plots:
        logger.info(f'Plotting top {top_n} annotations per position.')
        if svg:
            all_plot_path = ensure_path(output, f'{dataset_name}_cugo.svg', force=force)
        else:
            all_plot_path = ensure_path(output, f'{dataset_name}_cugo.png', force=force)

        # calculate dynamic figure width
        width = max(8, int((flank_upper - flank_lower + 1) * top_n * 0.6))
        figsize = (width, 16)

        # Create figure with gridspec for colorbar columns
        fig = plt.figure(figsize=figsize)
        gs = fig.add_gridspec(3, 2,
                              width_ratios=[20, 1],
                              height_ratios=[1, 1, 1],
                              hspace=0.7,
                              wspace=0.05)

        ax1 = fig.add_subplot(gs[0, 0])  # annotation plot
        ax2 = fig.add_subplot(gs[1, 0])  # Size plot
        cbar_ax1 = fig.add_subplot(gs[1, 1])  # Colorbar for size
        ax3 = fig.add_subplot(gs[2, 0])  # TMH plot
        cbar_ax2 = fig.add_subplot(gs[2, 1])  # Colorbar for TMH

        # plot top annotation per position
        pos_centers, pos_boundaries = plot_top_annotations_per_position(
            context_path=context_path,
            flank_lower=flank_lower,
            flank_upper=flank_upper,
            annotation=annotation,
            top_n=top_n,
            ax=ax1,
        )

        # plot protein size distribution
        norm_size, cmap_size = plot_size_per_position(
            context_path=context_path,
            flank_lower=flank_lower,
            flank_upper=flank_upper,
            title='Density of amino acid length per position',
            ax=ax2,
            pos_boundaries=pos_boundaries,
            bin_width=bin_width,
            y_range=y_range
        )

        # plot TMH distribution
        norm_tmh, cmap_tmh = plot_tmh_per_position(
            context_path=context_path,
            flank_lower=flank_lower,
            flank_upper=flank_upper,
            title='Density of transmembrane helix count per position',  # Fixed typo
            ax=ax3,
            pos_boundaries=pos_boundaries,
            y_range=tmh_y_range
        )

        ax1.tick_params(labelbottom=True)

        for ax in [ax1, ax2, ax3]:
            ax.set_aspect('auto')

        norm_fraction = Normalize(vmin=0, vmax=1)

        sm_size = ScalarMappable(norm=norm_fraction, cmap=cmap_size)
        sm_size.set_array([])
        cbar1 = plt.colorbar(sm_size, cax=cbar_ax1)
        cbar1.set_label('Relative density', fontsize=14)
        cbar1.ax.tick_params(labelsize=12)

        sm_tmh = ScalarMappable(norm=norm_fraction, cmap=cmap_tmh)
        sm_tmh.set_array([])
        cbar2 = plt.colorbar(sm_tmh, cax=cbar_ax2)
        cbar2.set_label('Relative density', fontsize=14)
        cbar2.ax.tick_params(labelsize=12)

        plt.savefig(all_plot_path, dpi=300, bbox_inches='tight')
        logger.info(f"Plot saved to {all_plot_path}")


# ======================================
# CUGO COMMAND LINE WORKFLOW
# ======================================

def cugo(db_path: str,
         cugo_range: int,
         fasta: str,
         id_list: str,
         output_dir: str,
         flank_lower: int,
         flank_upper: int,
         annotation: str = 'COG_ID',
         top_n: int = 3,
         threads: int = 1,
         svg: bool = False,
         force: bool = False,
         bin_width: int = 10,
         y_range: int = None,
         tmh_y_range: int = None):
    """
    Complete CUGO workflow: generate context data and create comprehensive plots.

    Args:
        db_path: Path to CUGO file
        cugo_range: Range around target protein for context
        output_dir: Directory for output files
        flank_lower: Lower flank boundary for plotting
        flank_upper: Upper flank boundary for plotting
        top_n: Number of top COGs to display (default: 3)
        threads: Number of threads (default: 1).
        svg: Generate plot in SVG format.
        force: Whether to overwrite existing files
        fasta: Optional path to FASTA file
        bin_width: Bin width for size plots
        y_range: Y-axis range for size plots

    Returns:
        tuple: (context_file_path, plot_file_path)
    """
    if annotation not in ANNOTATION_COLUMNS:
        raise ValueError(f'Invalid annotation. Please select one of the following annotations: {",".join(annotation_columns)}')

    # generate context data
    context_file = context(
        fasta=fasta,
        id_list=id_list,
        db_path=db_path,
        annotation=annotation,
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
        annotation=annotation,
        all_plots=True,
        bin_width=bin_width,
        y_range=y_range,
        tmh_y_range=tmh_y_range,
        svg=svg,
        force=force
    )

    # determine plot file path
    dataset_name = determine_dataset_name(context_file, '.', 0, '_context')
    if svg:
        plot_file = f"{dataset_name}_cugo.svg"
    else:
        plot_file = f"{dataset_name}_cugo.png"

    return context_file, plot_file

def cugo_select(context_path: str,
             position: int,
             db_path: str,
             output: str,
             force: bool = False):
    dataset_name = determine_dataset_name(context_path, '.', 0, '_context')
    output_path = ensure_path(output, f'{dataset_name}_{position}_ids.txt', force=force)

    # load context data
    df = pd.read_csv(context_path, sep='\t')

    # filter to specified position
    pos_data = df[df['position'] == position]

    # get unique seqIDs
    seq_ids = pos_data['seqID'].dropna().unique().tolist()

    if not seq_ids:
        logger.warning(f"No sequences found for position {position}")
        return output_path

    # connect to database
    conn = sqlite3.connect(db_path)

    # query database for sequences
    placeholders = ','.join('?' * len(seq_ids))
    cursor = conn.execute(f"""
        SELECT seqID, protein_seq 
        FROM protein_data 
        WHERE seqID IN ({placeholders}) AND protein_seq IS NOT NULL
    """, seq_ids)

    # write to file
    with open(output_path, 'w') as file:
        count = 0
        for seqid, compressed_seq in tqdm(cursor, total=len(seq_ids), desc="Retrieving sequences"):
            sequence = decompress_sequence(compressed_seq)
            file.write(f">{seqid}\n{sequence}\n")
            count += 1

    conn.close()
    logger.info(f"Retrieved {count} sequences to {output_path}")

    return output_path