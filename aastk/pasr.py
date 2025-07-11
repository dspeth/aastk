#!/usr/bin/env python3
from .util import determine_file_type, ensure_path, extract_unique_keys, read_fasta_to_dict, write_fa_matches, write_fq_matches

import subprocess
import matplotlib.pyplot as plt
import pandas as pd
import yaml
import logging
import sys
import json
from pathlib import Path
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


# default logger writes to log file needs implementing
logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.DEBUG,  # or logging.INFO if you want less verbosity
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    stream=sys.stdout
)
logging.getLogger('matplotlib').setLevel(logging.WARNING)
logging.getLogger("PIL").setLevel(logging.INFO)

def build_protein_db(protein_name: str,
                     seed_fasta: str,
                     threads: int,
                     db_dir: str,
                     force: bool = False):
    """
        Builds a DIAMOND protein database from a seed FASTA file.

        Args:
            protein_name (str): Name of the protein for the database.
            seed_fasta (str): Path to the FASTA file containing seed sequences.
            threads (int): Number of threads to use.
            db_dir (str): Directory where the database should be stored. (Default: current working directory)
            force (bool): If true, existing files/directories in output path are overwritten

        Returns:
            db_path: Path to DIAMOND protein database.

        Raises:
            RuntimeError: If the DIAMOND database creation fails.
    """
    # check for db_dir
    db_path = ensure_path(db_dir, f"{protein_name}_seed_db", force=force)

    # log the path
    logger.info(f"Building DIAMOND database for {protein_name} at {db_path}")

    # construct diamond command and output it for control
    cmd = ["diamond", "makedb",
           "--in", seed_fasta,
           "-d", db_path,
           "-p", str(threads)]
    logger.debug(f"Running command: {' '.join(cmd)}")

    #try running the subprocess for the Diamond makedb command
    try:
        result = subprocess.run(cmd, check=True)
        logger.debug(f"DIAMOND makedb output: {result.stdout}")

    except subprocess.CalledProcessError as e:
        logger.error(f"Error in building the DIAMOND database: {e}")
        logger.error(f"STDERR: {e.stderr}")
        raise RuntimeError(f"DIAMOND database creation failed: {e}") from e

    logger.info(f"Successfully built DIAMOND database at {db_path}")
    return db_path

def search_protein_db(db_path: str,
                      query_path: str,
                      protein_name: str,
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
        protein_name (str): Name of the target protein (used for output file naming or filtering).
        output_dir (str): Directory where results should be stored. (Default: current working directory)
        sensitivity (str): Choose sensitivity of diamond blastp search (Default: --fast)
        block (int): Choose diamond blastp sequence block size in billions of letters. (Default: 6)
        chunk (int): Choose number of chunks for diamond blastp index processing. (Default: 2)
        force (bool): If true, existing files/directories in output path are overwritten


    Returns:
        output_path: Path to tabular BLAST output file.
    """
    # check for output_dir
    output_path = ensure_path(output_dir, f"{protein_name}_hits.txt", force=force)
    column_info_path = ensure_path(output_dir, f"{protein_name}_columns.json", force=force)

    # define the output columns of interest
    columns = ["qseqid", "sseqid", "pident", "qlen", "slen", "length", "mismatch", "gapopen", "qstart", "qend",
               "sstart", "send", "evalue", "bitscore", "score"]

    # save column information to json file as to not hardcode the score column in the blast_score_ratio function
    column_info = {col: idx for idx, col in enumerate(columns)}
    with open(column_info_path, 'w') as f:
        json.dump(column_info, f, indent=2)

    logger.info(f"Saved column information to {column_info_path}")

    # check for sensitivity, if None set to default --fast
    sensitivity_param = f"--{sensitivity}" if sensitivity else "--fast"

    logger.info(f"Searching DIAMOND database {db_path} with query {query_path}")
    logger.info(f"Output path: {output_path}")
    logger.info(f"Using parameters: sensitivity={sensitivity_param}, block={block}, chunk={chunk}")

    try:
        # create blastp db search command
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

        # run diamond blastp
        subprocess.run(cmd, check=True)


    except subprocess.CalledProcessError as e:

        logger.error(f"Error in DIAMOND blastp search: {e}")
        logger.error(f"STDERR: {e.stderr}")
        raise RuntimeError(f"DIAMOND blastp search failed: {e}") from e

    logger.info(f"Successfully completed DIAMOND search. Results at {output_path}")
    return output_path, column_info_path

def extract_matching_sequences(protein_name: str,
                               blast_tab: str,
                               query_path: str,
                               output_dir: str,
                               key_column: int = 0,
                               force: bool = False):
    """
    Extracts reads that have BLAST/DIAMOND hits against a custom database.

    Args:
        protein_name (str): Name of the target protein (used for output file naming or filtering).
        blast_tab: Tabular BLAST/DIAMOND output file.
        query_path: Fasta or fastq file containing sequencing reads used as BLAST/DIAMOND queries.
        output_dir (str): Directory where extracted sequences should be stored.
        key_column: Column index in the BLAST tab file to pull unique IDs from (default is 0).
        force (bool): If true, existing files/directories in output path are overwritten

    Returns:
        out_fasta: Path to output FASTA file.
    """
    # check for output_dir
    out_fasta = ensure_path(output_dir, f"{protein_name}_matched.fasta", force=force)
    stats_path = ensure_path(output_dir, f"{protein_name}_matched.stats", force=force)

    # Extract unique keys (query IDs) from the specified column of the BLAST tab file
    matching_ids = extract_unique_keys(blast_tab, key_column)
    logger.info(f"Extracted {len(matching_ids)} unique matching sequence IDs")

    # Determine file type (fasta or fastq)
    file_type = determine_file_type(query_path)
    logger.info(f"Determined input file type: {file_type}")

    # keep track of number of sequences written to the output file; set up zero count
    sequences_written = 0

    # Open the output file for writing
    with open(out_fasta, "w") as out:
        # Use the appropriate generator based on the file type
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

    # also add seqkit stats output file to target dir
    cmd = ["seqkit", "stats",
            out_fasta,
            "-o", stats_path,
            "-a"]

    logger.debug(f"Running command: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

    return out_fasta, stats_path

def calculate_max_scores(protein_name: str,
                         extracted: str,
                         matrix: str,
                         output_dir: str,
                         force: bool = False):
    """
    Calculates max scores for sequences using a BLOSUM matrix.

    Args:
        protein_name (str): Name of the target protein (used for output file naming or filtering).
        extracted (str): Path to the extracted FASTA file.
        matrix (str): BLOSUM matrix name ('BLOSUM45' or 'BLOSUM62').
        output_dir (str): Directory where the output file should be stored.
        force (bool): If true, existing files/directories in output path are overwritten

    Returns:
        max_scores: Dictionary of protein headers and their max scores.
    """
    # check for output_dir
    out_file = ensure_path(output_dir, f"{protein_name}_max_scores.tsv", force=force)

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

    # read the fasta file containing matched sequences as obtained via pasr extract
    sequences = read_fasta_to_dict(extracted)
    logger.info(f"Read {len(sequences)} sequences from {extracted}")

    # calculate the max scores
    max_scores = {}
    for header, sequence in sequences.items():
        score = 0
        for amino_acid in sequence:
            aa = amino_acid.upper() # just make sure that we have consistency
            if aa in blosum_diagonals[matrix]:
                score += blosum_diagonals[matrix][amino_acid]
            else:
                logger.warning(f"Unknown amino acid '{amino_acid}' in sequence {header}")
        max_scores[header] = score

    # write the results of the step to specified output file
    with open(out_file, 'w') as out:
        out.write("Protein_id\tmax_score\n")
        for header, score in max_scores.items():
            out.write(f"{header}\t{score}\n")

    logger.info(f"Successfully wrote max scores to {out_file}")
    return out_file

def blast_score_ratio(protein_name: str,
                      blast_tab: str,
                      max_scores_path: str,
                      output_dir: str,
                      key_column: int = 0,
                      column_info_path: str = None,
                      score_column: int = None,
                      force: bool = False):
    """
    Computes BSR (Blast Score Ratio) using a BLAST tab file and max scores from a TSV.

    Args:
        protein_name (str): Name of the protein of interest.
        blast_tab (str): Path to DIAMOND/BLAST output file (must include 'score' column).
        max_scores_path (str): Path to TSV file with max scores (headers: Protein_id, max_score).
        output_dir (str): Directory to save the BSR results.
        key_column (int): Column index in blast_tab to use for matching. (Default: 0 for 'qseqid')
        force (bool): If true, existing files/directories in output path are overwritten

    Returns:
        bsr_output (str): Path to the output file with BSR values.
    """
    # check for output_dir
    bsr_file = ensure_path(output_dir, f"{protein_name}_bsr.tsv", force=force)

    logger.info(f"Computing blast score ratio (BSR) for {protein_name}")
    logger.info(f"Using blast tab file: {blast_tab}")
    logger.info(f"Using max scores file: {max_scores_path}")

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
    elif score_column:
        raw_score_column = score_column - 1
        logger.info(f"Using score column index {raw_score_column} from column info file")
    else:
        logger.warning(f"Column info file not provided or not found. "
                       f"Please enter path to column info file or pass "
                       f"the index of the score column in your BLAST tabular output")

    # parse the max_scores.tsv file
    max_scores = {}
    with open(max_scores_path) as tsv:
        header = tsv.readline()
        for line in tsv:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                max_scores[parts[0]] = float(parts[1])

    logger.info(f"Loaded {len(max_scores)} max scores")

    # process blast tab file and calculate BSR
    processed_count = 0
    error_count = 0

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
                raw_score = float(parts[raw_score_column])
                pident = float(parts[pident_column])
                qlen = int(parts[qlen_column])

                if key not in max_scores:
                    logger.warning(f"Line {line_num}: No max score found for '{key}'")
                    error_count += 1

                max_score = max_scores[key]
                if max_score <= 0:
                    logger.warning(f"Line {line_num}: Max score is zero or negative for '{key}': {max_score}")
                    error_count += 1
                    continue

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

def plot_bsr(protein_name: str,
             bsr_file: str,
             output_dir: str,
             yaml_path: str,
             force: bool = False,
             update: bool = False):
    """
    Creates a scatter plot of the BSR data.

    Args:
        protein_name: Name of the protein of interest.
        bsr_file: Path to the BSR TSV file.
        output_dir: Directory to save the plot.
        force (bool): If true, existing files/directories in output path are overwritten

    Returns:
        Path to the output plot.
    """
    # Determine output path
    if update:
        out_graph = ensure_path(output_dir, f'{protein_name}_updated_bsr.png', force=force)
    else:
        out_graph = ensure_path(output_dir, f'{protein_name}_bsr.png', force=force)

    logger.info(f"Creating BSR scatter plot for {protein_name}")

    try:
        bsr_df = pd.read_csv(bsr_file, sep='\t', header=0)

        # Check required columns
        required_cols = ['max_score', 'score', 'pident', 'qlen', 'BSR']
        for col in required_cols:
            if col not in bsr_df.columns:
                raise ValueError(f"Required column '{col}' not found in BSR data")

        bin_width = 10
        max_score_max = bsr_df['max_score'].max()
        score_max = bsr_df['score'].max()

        max_score_bins = np.arange(0, max_score_max + bin_width, bin_width)
        score_bins = np.arange(0, score_max + bin_width, bin_width)

        # Bin values
        bsr_df['max_score_bin'] = pd.cut(bsr_df['max_score'], bins=max_score_bins, include_lowest=True, right=False)
        bsr_df['score_bin'] = pd.cut(bsr_df['score'], bins=score_bins, include_lowest=True, right=False)

        # Count occurrences in each bin
        binned_counts = bsr_df.groupby(['max_score_bin', 'score_bin'], observed=True).size().reset_index(name='counts')

        # Extract midpoints of bins for positioning
        def bin_mid(bin_series):
            return bin_series.apply(lambda b: (b.left + b.right) / 2)

        binned_counts['x'] = bin_mid(binned_counts['max_score_bin'])
        binned_counts['y'] = bin_mid(binned_counts['score_bin'])

        # Create layout
        fig, axs = plt.subplot_mosaic(
            [['histx', '.'],
             ['scatter', 'histy']],
            figsize=(8, 8),
            width_ratios=(4, 1),
            height_ratios=(1, 4),
            layout='constrained'
        )

        # Scatterplot
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
        axs['scatter'].set_xlim(0, 1.5 * score_max)
        axs['scatter'].set_ylim(bottom=0)

        # Colorbar inside scatterplot
        cb_ax = inset_axes(axs['scatter'], width="5%", height="30%", loc='upper left', borderpad=1)
        cbar = fig.colorbar(scatter, cax=cb_ax)
        cbar.set_label('% sequence identity')

        # Histogram data using binned_counts with proper axis alignment
        # Top histogram (histx) - sum counts for each x position, align x-axis with scatter
        x_hist = binned_counts.groupby('x', observed=True)['counts'].sum()
        axs['histx'].bar(x_hist.index, x_hist.values, width=bin_width, align='center', color='gray', edgecolor='black')
        axs['histx'].set_xlim(axs['scatter'].get_xlim())
        axs['histx'].set_ylabel('Counts')
        axs['histx'].tick_params(labelbottom=False)
        axs['histx'].set_title(f'Protein Alignment Score Ratio for {protein_name}')

        y_hist = binned_counts.groupby('y', observed=True)['counts'].sum()
        axs['histy'].barh(y_hist.index, y_hist.values, height=bin_width, align='center', color='gray', edgecolor='black')
        axs['histy'].set_ylim(axs['scatter'].get_ylim())
        axs['histy'].set_xlabel('Counts')
        axs['histy'].tick_params(labelleft=False)

        max_count = y_hist.values.max() if len(y_hist) > 0 else 100
        axs['histy'].set_xticks(range(0, int(max_count) + 50, 50))

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
             protein_name: str = None,
             force: bool = False):
    """
    Writes metadata parameters to a YAML file.

    Args:
        selfmin: Minimum self-score threshold.
        selfmax: Maximum self-score threshold.
        output_dir: Directory to save the metadata file.
        dbmin: Minimum database score threshold.
        bsr: Minimum BSR threshold.
        protein_name (str): Name of the protein of interest.
        force (bool): If true, existing files/directories in output path are overwritten

    Returns:
        Path to the metadata file.
    """
    yaml_path = ensure_path(output_dir, f"{protein_name}.yaml", force=force)
    logger.info(f"Writing metadata for {protein_name}")

    params = {
        "protein_name": protein_name,
        "selfmin": selfmin,
        "selfmax": selfmax,
        "dbmin": dbmin,
        "bsr": bsr
    }

    # remove None values
    params = {k: v for k, v in params.items() if v is not None}
    yaml_str = yaml.dump(params)

    with open(yaml_path, "w") as f:
        f.write(yaml_str)

    logger.info(f"Successfully wrote metadata to {yaml_path}")
    return yaml_path


def select(yaml_path: str,
           matched_fasta: str,
           bsr_table: str,
           output_dir: str,
           selfmin: int = None,
           selfmax: int = None,
           dbmin: int = None,
           bsr: float = None,
           protein_name: str = None,
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
        protein_name (str): Name of the protein (required when params=True).
        force (bool): If true, existing files/directories in output path are overwritten.
        create_yaml (bool): Create a YAML file with the provided parameters.
        params (bool): Use provided parameters instead of YAML file (mutually exclusive with yaml_path).
    """
    logger.info("Starting sequence subsetting based on thresholds")

    # check mutual exclusivity of yaml_path and params
    if yaml_path is not None and params:
        raise ValueError(
            "yaml_path and params are mutually exclusive. Use either a YAML file or provide parameters directly.")

    if yaml_path is None and not params:
        raise ValueError("Either yaml_path must be provided or params must be True.")

    created_yaml_path = None

    if params:
        # validate required parameters when using params=True
        if selfmin is None or selfmax is None:
            raise ValueError("selfmin and selfmax are required when params=True")

        if protein_name is None:
            raise ValueError("protein_name is required when params=True")

        # ensure either dbmin or bsr is provided (but not both)
        if dbmin is None and bsr is None:
            raise ValueError("Either dbmin or bsr must be provided when params=True")

        if dbmin is not None and bsr is not None:
            raise ValueError("dbmin and bsr are mutually exclusive. Provide only one.")

        logger.info("Using provided parameters for thresholds")

        # Set the threshold values
        bsr_min = bsr

        # Create YAML file if requested
        if create_yaml:
            logger.info("Creating YAML file with provided parameters")
            created_yaml_path = metadata(
                selfmin=selfmin,
                selfmax=selfmax,
                output_dir=output_dir,
                dbmin=dbmin,
                bsr=bsr,
                protein_name=protein_name,
                force=force
            )

    else:
        # Load thresholds from YAML file
        logger.info(f"Loading thresholds from YAML file: {yaml_path}")
        try:
            with open(yaml_path) as f:
                thresholds = yaml.safe_load(f)
        except Exception as e:
            logger.error(f"Failed to load YAML file {yaml_path}: {e}")
            raise RuntimeError(f"Failed to load thresholds from {yaml_path}: {e}") from e

        # Get thresholds and assign default values
        protein_name = thresholds.get("protein_name", None)
        selfmin = thresholds.get("selfmin", 0)
        selfmax = thresholds.get("selfmax", float('inf'))
        dbmin = thresholds.get("dbmin", None)
        bsr_min = thresholds.get("bsr", None)

        if protein_name is None:
            raise ValueError("protein_name must be specified in the YAML file")

    # determine output path
    output_file = f"{protein_name}_matched_update.faa"
    stats_file = f"{protein_name}_matched_update.stats"

    output_path = ensure_path(output_dir, output_file, force=force)
    stats_path = ensure_path(output_dir, stats_file, force=force)

    # load BSR table
    bsr_df = pd.read_csv(bsr_table, sep='\t')

    # apply essential filters
    filtered = bsr_df[(bsr_df['max_score'] >= selfmin) & (bsr_df['max_score'] <= selfmax)]

    # apply mutually exclusive filters
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

    # load matched FASTA
    sequences = read_fasta_to_dict(matched_fasta)

    # track how many sequences we are keeping
    kept = 0

    # write subset
    with open(output_path, 'w') as out:
        for header, seq in sequences.items():
            if header in filtered_ids:
                out.write(f">{header}\n{seq}\n")
                kept += 1

    logger.info(f"Filtered sequences written to {output_path}: "
                f"{kept} sequences out of {len(sequences)} original sequences")

    # generate seqkit stats
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

def pasr(protein_name: str,
         seed_fasta: str,
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
        protein_name (str): Protein name.
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
    if update and not yaml_path:
        logger.error("YAML path is required if update is True")
        exit()

    # check for output_dir
    output_path = ensure_path(output_dir, force=force)
    logger.info(f"Running PASR workflow for {protein_name}")
    logger.info(f"Output directory: {output_path}")

    # store all the output paths in a dictionary
    results = {}

    try:
        logger.info("Building protein database")
        db_path = build_protein_db(protein_name, seed_fasta, threads, output_dir, force=force)
        results['db_path'] = db_path

        logger.info("Searching protein database")
        search_output, column_info_path = search_protein_db(db_path, query_fasta, protein_name, threads, output_path, sensitivity, block, chunk, force=force)
        results['search_output'] = search_output
        results['column_info_path'] = column_info_path

        logger.info("Extracting matching sequences")
        matched_fasta, stats_path = extract_matching_sequences(protein_name, search_output, query_fasta, output_path, key_column, force=force)
        results['matched_fasta'] = matched_fasta
        results['stats_path'] = stats_path

        logger.info("Calculating max scores")
        max_scores = calculate_max_scores(protein_name, matched_fasta, matrix_name, output_path, force=force)
        results['max_scores'] = max_scores

        logger.info("Calculating blast score ratios")
        bsr_file = blast_score_ratio(protein_name, search_output, max_scores, output_path, key_column, column_info_path, score_column=None, force=force)
        results['bsr_file'] = bsr_file

        logger.info("Creating BSR plot")
        bsr_plot = plot_bsr(protein_name, bsr_file, output_path, yaml_path, force=force, update=False)
        results['bsr_plot'] = bsr_plot

        if update:
            logger.info("Running update for specified data")
            subset_fasta, update_stats_path = select(yaml_path, matched_fasta, bsr_file, output_dir, force=force)
            results['subset_fasta'] = subset_fasta
            results['update_stats_path'] = update_stats_path
            updated_plot = plot_bsr(protein_name, bsr_file, output_path,  yaml_path, force=force, update=update)
            results['updated_plot'] = updated_plot

        logger.info("PASR workflow completed successfully")
        return results

    except Exception as e:
        logger.error(f"PASR workflow failed: {e}")
        exit()
