#!/usr/bin/env python3
from .util import *

import subprocess
import matplotlib.pyplot as plt
import pandas as pd
import yaml
import logging
import json
from pathlib import Path
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

logger = logging.getLogger(__name__)

def build_protein_db(seed_fasta: str,
                     threads: int,
                     db_dir: str,
                     force: bool = False):
    """
    Builds a DIAMOND protein database from a seed FASTA file.

    Args:
        seed_fasta (str): Path to the FASTA file containing seed sequences.
        threads (int): Number of threads to use.
        db_dir (str): Directory where the database should be stored. (Default: current working directory)
        force (bool): If true, existing files/directories in output path are overwritten

    Returns:
        db_path: Path to DIAMOND protein database.

    Raises:
        RuntimeError: If the DIAMOND database creation fails.
    """
    # check if diamond is in path
    check_dependency_availability('diamond')

    seed_fasta_filename = Path(seed_fasta).name
    protein_name = determine_dataset_name(seed_fasta_filename, '.', 0)

    # ===============================
    # Database path setup
    # ===============================
    db_path = ensure_path(db_dir, f"{protein_name}_seed_db", force=force)

    # log the path
    logger.info(f"Building DIAMOND database for {protein_name} at {db_path}")

    # =======================================
    # DIAMOND makedb command construction
    # =======================================
    cmd = ["diamond", "makedb",
           "--in", seed_fasta,
           "-d", db_path,
           "-p", str(threads)]
    logger.info(f"Running command: {' '.join(cmd)}")

    # ===============================
    # Execute DIAMOND makedb
    # ===============================
    try:
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            text=True
        )
        _, stderr = proc.communicate()

        if stderr:
            logger.log(99, stderr)

    except subprocess.CalledProcessError as e:
        logger.error(f"Error in building the DIAMOND database: {e}")
        logger.error(f"STDERR: {e.stderr}")
        raise RuntimeError(f"DIAMOND database creation failed: {e}") from e

    logger.info(f"Successfully built DIAMOND database at {db_path}")
    return db_path

def search_protein_db(db_path: str,
                      query_path: str,
                      threads: int,
                      output_dir: str,
                      sensitivity: str,
                      block: int = 6,
                      chunk: int = 2,
                      force: bool = False):
    """
    Searches a DIAMOND reference database for homologous sequences.

    Args:
        db_path (str): Path to the DIAMOND reference database.
        query_path (str): Path to the query FASTA file containing sequences to be searched.
        threads (int): Number of CPU threads to use for the search.
        output_dir (str): Directory where results should be stored. (Default: current working directory)
        sensitivity (str): Choose sensitivity of diamond blastp search (Default: --fast)
        block (int): Choose diamond blastp sequence block size in billions of letters. (Default: 6)
        chunk (int): Choose number of chunks for diamond blastp index processing. (Default: 2)
        force (bool): If true, existing files/directories in output path are overwritten


    Returns:
        output_path: Path to tabular BLAST output file.
    """
    # check if diamond is in PATH
    check_dependency_availability('diamond')

    # automatically name files
    db_name = Path(db_path).name
    protein_name = db_name.replace('_seed_db.dmnd', '')

    # ===============================
    # Output file path setup
    # ===============================
    output_path = ensure_path(output_dir, f"{protein_name}_hits.txt", force=force)
    column_info_path = ensure_path(output_dir, f"{protein_name}_columns.json", force=force)

    # =======================================================
    # DIAMOND blastp output file column configuration
    # =======================================================
    # define the output columns of interest
    columns = ["qseqid", "sseqid", "pident", "qlen", "slen", "length", "mismatch", "gapopen", "qstart", "qend",
               "sstart", "send", "evalue", "bitscore", "score"]

    # save column information to json file as to not hardcode the score column in the blast_score_ratio function
    column_info = {col: idx for idx, col in enumerate(columns)}
    with open(column_info_path, 'w') as f:
        json.dump(column_info, f, indent=2)

    logger.info(f"Saved column information to {column_info_path}")

    # ===============================
    # Parameter setup
    # ===============================
    # check for sensitivity, if None set to default --fast
    sensitivity_param = f"--{sensitivity}" if sensitivity else "--fast"

    logger.info(f"Searching DIAMOND database {db_path} with query {query_path}")
    logger.info(f"Output path: {output_path}")
    logger.info(f"Using parameters: sensitivity={sensitivity_param}, block={block}, chunk={chunk}")

    # =======================================
    # DIAMOND blastp command construction
    # =======================================
    cmd = ["diamond", "blastp",
           "-d", db_path,
           "-q", query_path,
           "-p", str(threads),
           "-o", output_path,
           "-k", str(1),
           "--matrix", "blosum45",
           "--masking", str(0),
           "--outfmt", str(6), *columns,
           "-b", str(block),
           "-c", str(chunk),
           "--min-score", str(50),
           "--comp-based-stats", str(0),
           sensitivity_param]

    logger.debug(f"Running command: {' '.join(cmd)}")
    # =======================================
    # Execute DIAMOND blastp search
    # =======================================
    try:
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            text=True
        )
        _, stderr = proc.communicate()

        if stderr:
            logger.log(99, stderr)
    except subprocess.CalledProcessError as e:

        logger.error(f"Error in DIAMOND blastp search: {e}")
        logger.error(f"STDERR: {e.stderr}")
        raise RuntimeError(f"DIAMOND blastp search failed: {e}") from e

    logger.info(f"Successfully completed DIAMOND search. Results at {output_path}")
    return output_path, column_info_path

def extract_matching_sequences(blast_tab: str,
                               query_path: str,
                               output_dir: str,
                               key_column: int = 0,
                               force: bool = False):
    """
    Extracts reads that have BLAST/DIAMOND hits against a custom database.

    Args:
        blast_tab: Tabular BLAST/DIAMOND output file.
        query_path: Fasta or fastq file containing sequencing reads used as BLAST/DIAMOND queries.
        output_dir (str): Directory where extracted sequences should be stored.
        key_column: Column index in the BLAST tab file to pull unique IDs from (default is 0).
        force (bool): If true, existing files/directories in output path are overwritten

    Returns:
        out_fasta: Path to output FASTA file.
    """
    # check if seqkit is in path
    check_dependency_availability('seqkit')

    blast_file = Path(blast_tab).name
    protein_name = blast_file.replace('_hits.txt', '')

    # ===============================
    # Output file setup
    # ===============================
    out_fasta = ensure_path(output_dir, f"{protein_name}_matched.faa", force=force)
    stats_path = ensure_path(output_dir, f"{protein_name}_matched.stats", force=force)

    # =========================================
    # Extract matching IDs from BLAST results
    # =========================================
    matching_ids = extract_unique_keys(blast_tab, key_column)
    logger.info(f"Extracted {len(matching_ids)} unique matching sequence IDs")

    # ===============================
    # Input file type detection
    # ===============================
    file_type = determine_file_type(query_path)
    logger.info(f"Determined input file type: {file_type}")

    # ======================================
    # Extract and write matching sequences
    # ======================================
    # keep track of number of sequences written to the output file; set up zero count
    sequences_written = 0

    with open(out_fasta, "w") as out:
        # use the appropriate generator based on the file type
        if file_type == "fasta":
            logger.info(f"Processing FASTA format from {query_path}")
            for header, sequence in write_fa_matches(query_path, matching_ids):
                out.write(f"{header}\n{sequence}\n")
                sequences_written += 1
        elif file_type == "fastq":
            logger.info(f"Processing FASTQ format from {query_path}")
            for header, sequence in write_fq_matches(query_path, matching_ids):
                out.write(f"{header}\n{sequence}\n")
                sequences_written += 1
        else:
            logger.error(f"Unsupported file type: {file_type}")
            raise ValueError(f"Unsupported file type: {file_type}")

        logger.info(f"Successfully wrote {sequences_written} matching sequences to {out_fasta}")

    # ===============================
    # Generate sequence statistics
    # ===============================
    cmd = ["seqkit", "stats",
            out_fasta,
            "-o", stats_path,
            "-a"]

    logger.debug(f"Running command: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

    return out_fasta, stats_path


def calculate_max_scores(extracted: str,
                         matrix: str,
                         output_dir: str,
                         force: bool = False):
    """
    Calculates max scores for sequences using a BLOSUM matrix.

    Args:
        extracted (str): Path to the extracted FASTA file.
        matrix (str): BLOSUM matrix name ('BLOSUM45' or 'BLOSUM62').
        output_dir (str): Directory where the output file should be stored.
        force (bool): If true, existing files/directories in output path are overwritten

    Returns:
        max_scores: Dictionary of protein headers and their max scores.
    """
    extracted_name = Path(extracted).name
    protein_name = extracted_name.replace('_matched.faa', '')

    # ===============================
    # Output file setup
    # ===============================
    out_file = ensure_path(output_dir, f"{protein_name}_max_scores.tsv", force=force)

    # ============================================================
    # BLOSUM matrix choice validation and BLOSUM diagonals setup
    # ============================================================
    valid_matrices = ["BLOSUM45", "BLOSUM62"]

    # check if valid matrix was chosen
    if matrix not in valid_matrices:
        logger.error(f"Invalid matrix: {matrix}. Must be one of {valid_matrices}")

    # define the diagonals of the usable BLOSUM matrices
    blosum_diagonals = {
        "BLOSUM45": {
            'A': 5, 'R': 7, 'N': 6, 'D': 7, 'C': 12, 'Q': 6, 'E': 6, 'G': 7,
            'H': 10, 'I': 5, 'L': 5, 'K': 5, 'M': 6, 'F': 8, 'P': 9, 'S': 4,
            'T': 5, 'W': 15, 'Y': 8, 'V': 5, 'B': 5, 'J': 4, 'Z': 5, 'X': 0
        },
        "BLOSUM62": {
            'A': 4, 'R': 5, 'N': 6, 'D': 6, 'C': 9, 'Q': 5, 'E': 5, 'G': 6,
            'H': 8, 'I': 4, 'L': 4, 'K': 5, 'M': 5, 'F': 6, 'P': 7, 'S': 4,
            'T': 5, 'W': 11, 'Y': 7, 'V': 4, 'B': 4, 'J': 3, 'Z': 4, 'X': 0
        }
    }

    logger.info(f"Calculating max scores using {matrix} matrix")

    # ===============================
    # Read input sequences
    # ===============================
    sequences = read_fasta_to_dict(extracted)
    logger.info(f"Read {len(sequences)} sequences from {extracted}")

    # ===============================
    # Calculate max scores
    # ===============================
    max_scores = {}
    for header, sequence in sequences.items():
        score = 0
        for amino_acid in sequence:
            aa = amino_acid.upper()  # just to make sure that we have consistency
            if aa in blosum_diagonals[matrix]:
                score += blosum_diagonals[matrix][amino_acid]
            else:
                logger.warning(f"Unknown amino acid '{amino_acid}' in sequence {header}")
        max_scores[header] = score

    # ===============================
    # Write results to output file
    # ===============================
    with open(out_file, 'w') as out:
        out.write("Protein_id\tmax_score\n")
        for header, score in max_scores.items():
            out.write(f"{header}\t{score}\n")

    logger.info(f"Successfully wrote max scores to {out_file}")

    return out_file

def blast_score_ratio(blast_tab: str,
                      max_scores_path: str,
                      output_dir: str,
                      key_column: int = 0,
                      column_info_path: str = None,
                      score_column: int = None,
                      force: bool = False):
    """
    Computes BSR (Blast Score Ratio) using a BLAST tab file and max scores from a TSV.

    Args:
        blast_tab (str): Path to DIAMOND/BLAST output file (must include 'score' column).
        max_scores_path (str): Path to TSV file with max scores (headers: Protein_id, max_score).
        output_dir (str): Directory to save the BSR results.
        key_column (int): Column index in blast_tab to use for matching. (Default: 0 for 'qseqid')
        column_info_path (str): Path to file containing column names for the blast tabular output file to retrieve index of the score column
        score_column (int): Index of blast tabular output file's score column
        force (bool): If true, existing files/directories in output path are overwritten

    Returns:
        bsr_output (str): Path to the output file with BSR values.
    """
    blast_file = Path(blast_tab).name
    protein_name = blast_file.replace('_hits.txt', '')

    # ===============================
    # Output file setup
    # ===============================
    bsr_file = ensure_path(output_dir, f"{protein_name}_bsr.tsv", force=force)

    logger.info(f"Computing blast score ratio (BSR) for {protein_name}")
    logger.info(f"Using blast tab file: {blast_tab}")
    logger.info(f"Using max scores file: {max_scores_path}")

    # ===============================================================================================
    # Determination of BLAST output column indices to be used for retrieval of relevant parameters
    # ===============================================================================================
    # retrieve indices for essential columns from column info file if path is provided
    if column_info_path and Path(column_info_path).exists():
        try:
            with open(column_info_path, 'r') as f:
                column_info = json.load(f)
                if 'score' in column_info:
                    raw_score_column = column_info['score']
                    logger.info(f"Using score column index {raw_score_column} from column info file")
                else:
                    logger.warning(f"'Score column not found in column info")
                if 'qlen' in column_info:
                    qlen_column = column_info['qlen']
                    logger.info(f"Using qlen column index {qlen_column} from column info file")
                else:
                    logger.warning(f"'qlen column not found in column info")
                if 'pident' in column_info:
                    pident_column = column_info['pident']
                    logger.info(f"Using pident column index {pident_column} from column info file")
                else:
                    logger.warning(f"'pident column not found in column info")
        except (json.JSONDecodeError, IOError) as e:
            logger.warning(f"Error reading column info file: {e}")
    # if no column info file is provided we use the input score column index
    elif score_column:
        raw_score_column = score_column - 1
        logger.info(f"Using score column index {raw_score_column} from column info file")
    else:
        logger.warning(f"Column info file not provided or not found. "
                       f"Please enter path to column info file or pass "
                       f"the index of the score column in your BLAST tabular output")

    # ===============================
    # Load max scores data
    # ===============================
    max_scores = {}
    with open(max_scores_path) as tsv:
        header = tsv.readline()
        for line in tsv:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                max_scores[parts[0]] = float(parts[1])

    logger.info(f"Loaded {len(max_scores)} max scores")

    # ===============================
    # Process BLAST results and calculate BSR
    # ===============================
    processed_count = 0
    error_count = 0

    # write the blast score ratio output file
    with open(blast_tab) as infile, open(bsr_file, 'w') as out:
        out.write("qseqid\tsseqid\tpident\tqlen\tscore\tmax_score\tBSR\n")

        for line_num, line in enumerate(infile, 1):
            parts = line.strip().split('\t')

            if len(parts) <= max(key_column, raw_score_column):
                logger.warning(f"Line {line_num}: Not enough columns "
                               f"(expected at least {max(key_column, score_column)+1}, got {len(parts)})")
                error_count += 1
                continue

            key = parts[key_column]
            try:
                # raw score is the alignment score
                raw_score = float(parts[raw_score_column])
                # % sequence identity (used to color the datapoints in the plot later)
                pident = float(parts[pident_column])
                # query length as useful additional information
                qlen = int(parts[qlen_column])

                if key not in max_scores:
                    logger.warning(f"Line {line_num}: No max score found for '{key}'")
                    error_count += 1

                max_score = max_scores[key]
                if max_score <= 0:
                    logger.warning(f"Line {line_num}: Max score is zero or negative for '{key}': {max_score}")
                    error_count += 1
                    continue

                # calculate the blast score ratio
                bsr = raw_score / max_score
                out.write(f"{parts[0]}\t{parts[1]}\t{pident:.2f}\t{qlen}\t{raw_score:.1f}\t{max_score:.1f}\t{bsr:.4f}\n")
                processed_count += 1

            except (KeyError, ValueError, IndexError) as e:
                logger.warning(f"Line {line_num}: Error processing '{line.strip()}': {e}")
                if isinstance(e, IndexError):
                    logger.warning(
                        f"Line {line_num}: Tried to access column {score_column} but line only has {len(parts)} columns")
                error_count += 1

    logger.info(f"Successfully processed {processed_count} entries and wrote BSR to {bsr_file}")
    if error_count > 0:
        logger.warning(f"Encountered {error_count} errors during BSR calculation")
    return bsr_file

def plot_bsr(bsr_file: str,
             output_dir: str,
             yaml_path: str,
             force: bool = False,
             update: bool = False):
    """
    Creates a scatter plot of the BSR data flanked by histograms showing the distribution of datapoints alongside the axes.

    Args:
        bsr_file (str): Path to the BSR TSV file.
        output_dir (str): Directory to save the plot.
        yaml_path (str): Path to yaml file containing metadata which is plotted if update is selected
        force (bool): If true, existing files/directories in output path are overwritten
        update (bool): If true, metadata from the yaml path is plotted within the normal BSR plot.

    Returns:
        Path to the output plot.
    """
    bsr_file_name = Path(bsr_file).name
    protein_name = bsr_file_name.replace('_bsr.tsv', '')

    # ===============================
    # Output file path setup
    # ===============================
    if update:
        out_graph = ensure_path(output_dir, f'{protein_name}_updated_bsr.png', force=force)
    else:
        out_graph = ensure_path(output_dir, f'{protein_name}_bsr.png', force=force)

    logger.info(f"Creating BSR scatter plot for {protein_name}")

    # ===============================
    # Data loading and validation
    # ===============================
    try:
        bsr_df = pd.read_csv(bsr_file, sep='\t', header=0)

        # check BSR file consistency
        required_cols = ['max_score', 'score', 'pident', 'qlen', 'BSR']
        for col in required_cols:
            if col not in bsr_df.columns:
                raise ValueError(f"Required column '{col}' not found in BSR data")

        # ===============================
        # Histogram binning setup
        # ===============================
        # set bin width for flanking histograms
        bin_width = 10

        # determine number of bins according to maximum theoretical max score and alignment score in the data; arange the bins
        max_score_max = bsr_df['max_score'].max()
        score_max = bsr_df['score'].max()
        max_score_bins = np.arange(0, max_score_max + bin_width, bin_width)
        score_bins = np.arange(0, score_max + bin_width, bin_width)

        # bin theoretical max scores as well as alignment scores
        bsr_df['max_score_bin'] = pd.cut(bsr_df['max_score'], bins=max_score_bins, include_lowest=True, right=False)
        bsr_df['score_bin'] = pd.cut(bsr_df['score'], bins=score_bins, include_lowest=True, right=False)

        # count occurrences of datapoints in each bin
        binned_counts = bsr_df.groupby(['max_score_bin', 'score_bin'], observed=True).size().reset_index(name='counts')

        binned_counts['x'] = bin_mid(binned_counts['max_score_bin'])
        binned_counts['y'] = bin_mid(binned_counts['score_bin'])

        # ===============================
        # Plot layout and scatter plot
        # ===============================
        # create plot layout
        fig, axs = plt.subplot_mosaic(
            [['histx', '.'],
            ['scatter', 'histy']],
            figsize=(8, 8),
            width_ratios=(4, 1),
            height_ratios=(1, 4),
            layout='constrained'
        )

        # scatterplot
        scatter = axs['scatter'].scatter(
            bsr_df['max_score'],
            bsr_df['score'],
            c=bsr_df['pident'],
            cmap='viridis',
            edgecolor='k',
            alpha=0.5,
            vmin=0,
            vmax=100
        )
        axs['scatter'].set_xlabel('Calculated maximum score')
        axs['scatter'].set_ylabel('Alignment score to seed set')
        axs['scatter'].set_xlim(0, 1.5 * score_max) # set to 1.5x the maximum alignment score to capture relevant data in the scatterplot
        axs['scatter'].set_ylim(bottom=0)

        # move colorbar inside scatterplot and color it by sequence identity
        cb_ax = inset_axes(axs['scatter'], width="5%", height="30%", loc='upper left', borderpad=1)
        cbar = fig.colorbar(scatter, cax=cb_ax)
        cbar.set_label('% sequence identity')

        # ===============================
        # Histogram creation
        # ===============================
        # histogram data using binned_counts with proper axis alignment
        # top histogram (histx) - sum counts for each x position, align x-axis with scatter
        x_hist = binned_counts.groupby('x', observed=True)['counts'].sum()
        axs['histx'].bar(x_hist.index, x_hist.values, width=bin_width, align='center', color='gray', edgecolor='black')
        axs['histx'].set_xlim(axs['scatter'].get_xlim())
        axs['histx'].set_xticks([0, max_score_max / 2, max_score_max]) # maybe change to 1.5x
        axs['histx'].set_ylabel('Counts')
        axs['histx'].tick_params(labelbottom=False)
        axs['histx'].set_title(f'Protein Alignment Score Ratio for {protein_name}')

        # right histogram (histy) - sum counts for each y position, align x-axis with scatterplot's y-axis
        y_hist = binned_counts.groupby('y', observed=True)['counts'].sum()
        axs['histy'].barh(y_hist.index, y_hist.values, height=bin_width, align='center', color='gray', edgecolor='black')
        axs['histy'].set_ylim(axs['scatter'].get_ylim())
        axs['histy'].set_xlabel('Counts')
        axs['histy'].set_yticks([0, score_max / 2, score_max])
        axs['histy'].tick_params(labelleft=False)

        max_count = y_hist.values.max() if len(y_hist) > 0 else 100
        axs['histy'].set_xticks(range(0, int(max_count) + 50, 50))

        # ===============================
        # Threshold lines (update mode)
        # ===============================
        # plot the parameters from the metadata yaml file in case of update
        if update:
            try:
                with open(yaml_path) as f:
                    thresholds = yaml.safe_load(f)
            except Exception as e:
                raise RuntimeError(f"Failed to load thresholds from {yaml_path}: {e}") from e

            selfmin = thresholds.get("selfmin", 0)
            selfmax = thresholds.get("selfmax", float('inf'))
            dbmin = thresholds.get("dbmin", None)
            bsr_min = thresholds.get("bsr", None)

            if selfmin is not None:
                axs['scatter'].axvline(selfmin, color='black', linestyle='--', linewidth=1.0)
            if selfmax is not None:
                axs['scatter'].axvline(selfmax, color='black', linestyle='--', linewidth=1.0)
            if dbmin is not None:
                axs['scatter'].axhline(dbmin, color='black', linestyle='--', linewidth=1.0)
            if bsr_min is not None:
                x_min, x_max = axs['scatter'].get_xlim()
                x_vals = np.linspace(x_min, x_max, 500)
                y_vals = bsr_min * x_vals
                axs['scatter'].plot(x_vals, y_vals, color='black', linestyle='--', linewidth=1.0)

        # ===============================
        # Save plot and cleanup
        # ===============================
        fig.savefig(out_graph, dpi=300)
        plt.close(fig)

        logger.info(f"Successfully created BSR plot at {out_graph}")
        return out_graph

    except Exception as e:
        logger.error(f"Error creating BSR plot: {e}")
        raise RuntimeError(f"Failed to create BSR plot: {e}") from e


def metadata(selfmin: int,
             selfmax: int,
             output_dir: str,
             dbmin: int = None,
             bsr: float = None,
             seed: str = None,
             bsr_file: str = None,
             force: bool = False):
    """
    Writes metadata parameters to a YAML file.

    Args:
        selfmin: Minimum self-score threshold.
        selfmax: Maximum self-score threshold.
        output_dir: Directory to save the metadata file.
        dbmin: Minimum database score threshold.
        bsr: Minimum BSR threshold.
        seed (str): Path to seed FASTA.
        bsr_file (str): Path to BSR tab separated file.
        force (bool): If true, existing files/directories in output path are overwritten.

    Returns:
        Path to the metadata file.
    """
    logger.info("Starting metadata creation")

    # ===============================
    # Parameter validation
    # ===============================
    if seed is None and bsr_file is None:
        error_msg = "Either seed or bsr_file must be provided"
        logger.error(error_msg)
        raise ValueError(error_msg)

    if seed is not None and bsr_file is not None:
        error_msg = "seed and bsr_file are mutually exclusive. Provide only one."
        logger.error(error_msg)
        raise ValueError(error_msg)

    # ===============================
    # Protein name determination
    # ===============================
    if seed is not None:
        logger.info(f"Determining protein name from seed file: {seed}")
        seed_name = Path(seed).name
        protein_name = determine_dataset_name(seed_name, '.', 0)
        logger.info(f"Extracted protein name '{protein_name}' from seed file")
    elif bsr_file is not None:
        logger.info(f"Determining protein name from BSR file: {bsr_file}")
        bsr_name = Path(bsr_file).name
        protein_name = bsr_name.replace('_bsr.tsv', '')
        logger.info(f"Extracted protein name '{protein_name}' from BSR file")
    else:
        error_msg = "Unexpected state: both seed and bsr_file are None"
        logger.error(error_msg)
        raise RuntimeError(error_msg)

    # ===============================
    # Output file path setup
    # ===============================
    try:
        yaml_path = ensure_path(output_dir, f"{protein_name}.yaml", force=force)
        logger.info(f"Writing metadata for {protein_name} to {yaml_path}")
    except Exception as e:
        logger.error(f"Failed to create output path: {e}")
        raise

    # ===============================
    # Parameter collection and validation
    # ===============================
    logger.debug(f"Collecting parameters: selfmin={selfmin}, selfmax={selfmax}, "
                 f"dbmin={dbmin}, bsr={bsr}")

    # Validate threshold values
    if selfmin < 0 or selfmax < 0:
        error_msg = "selfmin and selfmax must be non-negative"
        logger.error(error_msg)
        raise ValueError(error_msg)

    if selfmin > selfmax:
        error_msg = f"selfmin ({selfmin}) cannot be greater than selfmax ({selfmax})"
        logger.error(error_msg)
        raise ValueError(error_msg)

    if bsr is not None and (bsr < 0 or bsr > 1):
        error_msg = f"BSR value must be between 0 and 1, got {bsr}"
        logger.error(error_msg)
        raise ValueError(error_msg)

    if dbmin is not None and dbmin < 0:
        error_msg = f"dbmin must be non-negative, got {dbmin}"
        logger.error(error_msg)
        raise ValueError(error_msg)

    params = {
        "protein_name": protein_name,
        "selfmin": selfmin,
        "selfmax": selfmax,
        "dbmin": dbmin,
        "bsr": bsr
    }

    # remove None values
    params = {k: v for k, v in params.items() if v is not None}
    logger.debug(f"Final parameters to write: {params}")

    # ===============================
    # Write YAML file
    # ===============================
    try:
        yaml_str = yaml.dump(params)
        logger.debug(f"Generated YAML content:\n{yaml_str}")

        with open(yaml_path, "w") as f:
            f.write(yaml_str)

        logger.info(f"Successfully wrote metadata to {yaml_path}")
        return yaml_path

    except Exception as e:
        logger.error(f"Failed to write metadata file {yaml_path}: {e}")
        raise RuntimeError(f"Failed to write metadata file: {e}") from e

def select(yaml_path: str,
           matched_fasta: str,
           bsr_table: str,
           output_dir: str,
           selfmin: int = None,
           selfmax: int = None,
           dbmin: int = None,
           bsr: float = None,
           force: bool = False,
           create_yaml: bool = False,
           params: bool = False):
    """
    Subsets matched sequences based on YAML thresholds or provided parameters.

    Args:
        yaml_path(str): Path to the metadata yaml file.
        matched_fasta (str): Path to the matched sequences FASTA.
        bsr_table (str): Path to the BSR table.
        output_dir (str): Directory to save the subsetted sequences.
        selfmin (int): Minimum self score threshold (required when params=True).
        selfmax (int): Maximum self score threshold (required when params=True).
        dbmin (int): Minimum database score threshold (mutually exclusive with bsr).
        bsr (float): Minimum BSR threshold (mutually exclusive with dbmin).
        force (bool): If true, existing files/directories in output path are overwritten.
        create_yaml (bool): Create a YAML file with the provided parameters.
        params (bool): Use provided parameters instead of YAML file (mutually exclusive with yaml_path).
    """
    # check if seqkit is in path
    check_dependency_availability('seqkit')

    logger.info("Starting sequence subsetting based on thresholds")

    bsr_name = Path(bsr_table).name
    protein_name = bsr_name.replace('_bsr.tsv', '')

    # ===============================
    # Parameter validation
    # ===============================
    # check mutual exclusivity of yaml_path and params
    if yaml_path is not None and params:
        raise ValueError(
            "yaml_path and params are mutually exclusive. Use either a YAML file or provide parameters directly.")

    if yaml_path is None and not params:
        raise ValueError("Either yaml_path must be provided or params must be True.")

    created_yaml_path = None

    # ===============================
    # Threshold loading/validation
    # ===============================
    if params:
        # validate required parameters when using params=True
        if selfmin is None or selfmax is None:
            raise ValueError("selfmin and selfmax are required when params=True")

        # ensure either dbmin or bsr is provided (but not both)
        if dbmin is None and bsr is None:
            raise ValueError("Either dbmin or bsr must be provided when params=True")

        if dbmin is not None and bsr is not None:
            raise ValueError("dbmin and bsr are mutually exclusive. Provide only one.")

        logger.info("Using provided parameters for thresholds")

        # set the threshold values
        bsr_min = bsr

        # create YAML file if requested
        if create_yaml:
            logger.info("Creating YAML file with provided parameters")
            created_yaml_path = metadata(
                selfmin=selfmin,
                selfmax=selfmax,
                output_dir=output_dir,
                dbmin=dbmin,
                bsr=bsr,
                bsr_file=bsr_table,
                force=force
            )

    else:
        # load thresholds from YAML file
        logger.info(f"Loading thresholds from YAML file: {yaml_path}")
        try:
            with open(yaml_path) as f:
                thresholds = yaml.safe_load(f)
        except Exception as e:
            logger.error(f"Failed to load YAML file {yaml_path}: {e}")
            raise RuntimeError(f"Failed to load thresholds from {yaml_path}: {e}") from e

        # get thresholds and assign default values
        protein_name = thresholds.get("protein_name", None)
        selfmin = thresholds.get("selfmin", 0)
        selfmax = thresholds.get("selfmax", float('inf'))
        dbmin = thresholds.get("dbmin", None)
        bsr_min = thresholds.get("bsr", None)

        if protein_name is None:
            raise ValueError("protein_name must be specified in the YAML file")

    # ===============================
    # Output file path setup
    # ===============================
    # create output file path
    output_file = f"{protein_name}_matched_update.faa"
    stats_file = f"{protein_name}_matched_update.stats"
    output_path = ensure_path(output_dir, output_file, force=force)
    stats_path = ensure_path(output_dir, stats_file, force=force)

    # ===============================
    # Data loading and filtering
    # ===============================
    # load BSR table
    bsr_df = pd.read_csv(bsr_table, sep='\t')

    # apply filters defined by the required arguments selfmin and selfmax to the BSR table
    filtered = bsr_df[(bsr_df['max_score'] >= selfmin) & (bsr_df['max_score'] <= selfmax)]

    # apply mutually exclusive filters to the BSR table (either min. alignment score or min. BSR)
    if dbmin is not None:
        filtered = filtered[filtered['score'] >= dbmin]
        logger.info(f"Applied dbmin filter: score >= {dbmin}")
    elif bsr_min is not None:
        filtered = filtered[filtered['BSR'] >= bsr_min]
        logger.info(f"Applied BSR filter: BSR >= {bsr_min}")

    if filtered.empty:
        logger.warning("No sequences passed the filtering criteria!")

    # extract sequence IDs kept after filtering
    filtered_ids = set(filtered['qseqid'])

    # ===============================
    # Sequence filtering and writing
    # ===============================
    # load matched FASTA
    sequences = read_fasta_to_dict(matched_fasta)

    # track how many sequences we are keeping
    kept = 0

    # retrieve filtered sequences from original homologous faa file
    with open(output_path, 'w') as out:
        for header, seq in sequences.items():
            if header in filtered_ids:
                out.write(f">{header}\n{seq}\n")
                kept += 1

    logger.info(f"Filtered sequences written to {output_path}: "
                f"{kept} sequences out of {len(sequences)} original sequences")

    # ===============================
    # Statistics generation
    # ===============================
    # generate seqkit stats for the filtered fasta file
    cmd = ["seqkit", "stats",
           output_path,
           "-o", stats_path,
           "-a"]

    logger.debug(f"Running command: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

    if created_yaml_path:
        return output_path, stats_path, created_yaml_path
    else:
        return output_path, stats_path

def pasr(seed_fasta: str,
         query_fasta: str,
         matrix_name: str,
         output_dir: str,
         sensitivity: str,
         block: int,
         chunk: int,
         key_column: int = 0,
         threads: int = 1,
         update: bool = False,
         yaml_path: str = None,
         force: bool = False):
    """
    PASR workflow with configurable output directory.

    Args:
        seed_fasta (str): Path to seed FASTA.
        query_fasta (str): Path to query FASTA.
        matrix_name (str): BLOSUM matrix ('BLOSUM45' or 'BLOSUM62').
        output_dir (str): Output directory (default: current directory).
        sensitivity (str): Choose sensitivity of diamond blastp search (Default: --fast)
        block (int): Choose diamond blastp sequence block size in billions of letters. (Default: 6)
        chunk (int): Choose number of chunks for diamond blastp index processing. (Default: 2)
        key_column: Column index in the BLAST tab file to pull unique IDs from (default is 0).
        threads (int): Number of threads (default: 1).
        update (bool): Creates updated matched_fasta and bsr plot using metadata yaml file
        yaml_path (str): Path to metadata yaml file
        force (bool): If true, existing files/directories in output path are overwritten
    """
    # ===============================
    # Parameter validation and setup
    # ===============================
    if update and not yaml_path:
        logger.error("YAML path is required if update is True")
        exit()

    seed_name = Path(seed_fasta).name
    protein_name = determine_dataset_name(seed_name, '.', 0)

    # check for output_dir
    output_path = ensure_path(output_dir, force=force)
    logger.info(f"Running PASR workflow for {protein_name}")
    logger.info(f"Output directory: {output_path}")

    # store all the output paths in a dictionary
    results = {}

    # ===============================
    # Database building
    # ===============================
    try:
        logger.info("Building protein database")
        db_path = build_protein_db(seed_fasta, threads, output_dir, force=force)
        results['db_path'] = db_path

        # ===============================
        # Database search
        # ===============================
        logger.info("Searching protein database")
        search_output, column_info_path = search_protein_db(db_path, query_fasta, threads, output_path, sensitivity, block, chunk, force=force)
        results['search_output'] = search_output
        results['column_info_path'] = column_info_path

        # ===============================
        # Sequence extraction
        # ===============================
        logger.info("Extracting matching sequences")
        matched_fasta, stats_path = extract_matching_sequences(search_output, query_fasta, output_path, key_column, force=force)
        results['matched_fasta'] = matched_fasta
        results['stats_path'] = stats_path

        # ===============================
        # Score calculations
        # ===============================
        logger.info("Calculating max scores")
        max_scores = calculate_max_scores(matched_fasta, matrix_name, output_path, force=force)
        results['max_scores'] = max_scores

        logger.info("Calculating blast score ratios")
        bsr_file = blast_score_ratio(search_output, max_scores, output_path, key_column, column_info_path, score_column=None, force=force)
        results['bsr_file'] = bsr_file

        # ===============================
        # Visualization
        # ===============================
        logger.info("Creating BSR plot")
        bsr_plot = plot_bsr(bsr_file, output_path, yaml_path, force=force, update=False)
        results['bsr_plot'] = bsr_plot

        # ===============================
        # Update workflow (optional)
        # ===============================
        if update:
            logger.info("Running update for specified data")
            subset_fasta, update_stats_path = select(yaml_path, matched_fasta, bsr_file, output_dir, force=force)
            results['subset_fasta'] = subset_fasta
            results['update_stats_path'] = update_stats_path
            updated_plot = plot_bsr(bsr_file, output_path,  yaml_path, force=force, update=update)
            results['updated_plot'] = updated_plot

        logger.info("PASR workflow completed successfully")
        return results

    except Exception as e:
        logger.error(f"PASR workflow failed: {e}")
        exit()