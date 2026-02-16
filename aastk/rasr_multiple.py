#!/usr/bin/env python3
import json
from .util import *
from .pasr import build as pasr_build
import subprocess
import logging
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


logger = logging.getLogger(__name__)


# ===============================
# Input normalization function
# ===============================
def list_inputs(gene_db_fasta, query_fastq):
    """
    Normalizes inputs to lists of files
    Accepts both individual files and directories

    Args:
        gene_db_fasta (str): Path to FASTA file or directory of FASTA files
        query_fastq (str): Path to FASTQ file or directory of FASTQ files
    
    Returns:
        db_files (list): List of paths to FASTA files
        query_files (list): List of paths to FASTQ files
    """
    # Normalize gene database inputs
    if Path(gene_db_fasta).is_dir(): 
        db_files = []
        for ext in ["*.fasta", "*.fa", "*.faa"]: 
            db_files.extend(Path(gene_db_fasta).glob(ext)) 
        db_files.sort()  # Sort files for consistent order
    else: 
        db_files = [Path(gene_db_fasta)]
    
    # Normalize query inputs 
    if Path(query_fastq).is_dir(): 
        query_files = []
        for ext in ["*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz"]:     
            query_files.extend(Path(query_fastq).glob(ext)) 
        query_files.sort()  
    else: 
        query_files = [Path(query_fastq)]
    
    logger.info(f"Found {len(query_files)} query file(s)")
    logger.info(f"Found {len(db_files)} database file(s)")

    return db_files, query_files




# ===============================
# Databases merging Function
# ===============================
def merge_dbs(db_files: list,
              output_dir: str,
              force: bool = False):
    """
    Merges multiple FASTA files into a single unique FASTA file.

    Args:
        db_files (list): List of paths to FASTA files
        output_dir (str): Directory to save merged FASTA file
        force (bool): Whether to overwrite existing files

    Returns:
        merged_fasta_path (str): Path to merged FASTA file.
        db_dict (dict): Dictionary mapping db names to their sequence IDs.
    """
    merged_fasta_path = ensure_path(output_dir, "merged_gene_db.fasta", force=force)
    
    if Path(merged_fasta_path).exists() and not force:
        logger.info(f"Merged database already exists at {merged_fasta_path}, skipping merge.")
        
        # Reconstruct db_dict by reading input files
        db_dict = {}
        for db_file in db_files:
            gene_name = db_file.stem
            
            with open(db_file, 'r') as db_fasta:
                for line in db_fasta:
                    line = line.strip()
                    if line.startswith('>'):
                        header = line[1:]  # Remove '>'
                        if gene_name not in db_dict:
                            db_dict[gene_name] = []
                        db_dict[gene_name].append(header)
        
        logger.info(f"Reconstructed db_dict from {len(db_files)} input file(s)")
        return merged_fasta_path, db_dict
    
    with open(merged_fasta_path, 'w') as outfile:
        db_dict = {}
        
        for db_file in db_files:
            gene_name = db_file.stem
            
            with open(db_file, 'r') as db_fasta:
                current_header = None
                current_seq = []

                for line in db_fasta:
                    line = line.strip()
                    if line.startswith('>'):
                        # Write previous sequence if exists
                        if current_header:
                            outfile.write(f">{current_header}\n")
                            outfile.write(''.join(current_seq) + '\n')
                            if gene_name not in db_dict:
                                db_dict[gene_name] = []
                            db_dict[gene_name].append(current_header)
                        
                        # Start new sequence
                        current_header = line[1:]  # Remove '>'
                        current_seq = []
                    else:
                        current_seq.append(line)
                
                # Write last sequence
                if current_header:
                    outfile.write(f">{current_header}\n")
                    outfile.write(''.join(current_seq) + '\n')
                    if gene_name not in db_dict:
                        db_dict[gene_name] = []
                    db_dict[gene_name].append(current_header)
    
    logger.info(f"Merged database created at {merged_fasta_path}")
    return str(merged_fasta_path), db_dict




# ===============================
# aastk search CLI FUNCTION
# ===============================
def search(gene_db_out_path: str,
            query_fastq: str,
            threads: int,
            output_dir: str,
            sensitivity: str,
            block: int = 6,
            chunk: int = 2,
            force: bool = False):
    """
    Searches a DIAMOND reference database for homologous sequences.
    Single database vs. single query

    Args:
        gene_db_out_path (str): Path to DIAMOND protein database for the gene of interest
        query_fastq (str): Path to sequencing read file, can be gzipped
        threads (int): Number of threads to use
        output_dir (str): Directory to save output files
        sensitivity (str): Sensitivity setting for DIAMOND search
        block (int): Block size parameter for DIAMOND search
        chunk (int): Chunk size parameter for DIAMOND search
        force (bool): Whether to overwrite existing files

    Returns:
        blast_out_path: Path to tabular BLAST output file.
    """
    
    # Check if diamond is available
    check_dependency_availability("diamond")

    # automatically name files
    db_filename = determine_dataset_name(gene_db_out_path, '.', 0, '_db')

    # ===============================
    # Output file path setup
    # ===============================
    output_path = ensure_path(output_dir, f"{db_filename}_hits.tsv", force=force)
    column_info_path = ensure_path(output_dir, f"{db_filename}_columns.json", force=force)

    # =======================================================
    # DIAMOND blastp output file column configuration
    # =======================================================
    # define the output columns of interest
    columns = ["qseqid", "sseqid", "pident", "qlen", "slen", "length", "mismatch", "gapopen", "qstart", "qend",
               "sstart", "send", "evalue", "bitscore", "score"]

    # save column information to json file as to not hardcode the score column in the bsr function
    column_info = {col: idx for idx, col in enumerate(columns)}
    with open(column_info_path, 'w') as f:
        json.dump(column_info, f, indent=2)

    logger.info(f"Saved column information to {column_info_path}")

    # ===============================
    # Parameter setup
    # ===============================
    # check for sensitivity, if None set to default --fast
    sensitivity_param = f"--{sensitivity}" if sensitivity else "--fast"

    logger.info(f"Searching DIAMOND database {gene_db_out_path} with query {query_fastq}")
    logger.info(f"Output path: {output_path}")
    logger.info(f"Using parameters: sensitivity={sensitivity_param}, block={block}, chunk={chunk}")

    # =======================================
    # DIAMOND blastx command construction
    # =======================================
    min_score = 10
    cmd = ["diamond", "blastx",
           "-d", gene_db_out_path,
           "-q", query_fastq,
           "-p", str(threads),
           "-o", output_path,
           "--max-target-seqs", str(1),
           "--matrix", "blosum45",
           "--masking", str(0),
           "--outfmt", str(6), *columns,
           "-b", str(block),
           "-c", str(chunk),
           "--min-score", str(min_score),
           "--comp-based-stats", str(0),
           sensitivity_param]

    logger.debug(f"Running command: {' '.join(cmd)}")

    # =======================================
    # Execute DIAMOND blastx search
    # =======================================
    try:
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            text=True
        )
        _, stderr = proc.communicate()

        if proc.returncode != 0:
            logger.error(f"DIAMOND blastx failed with return code {proc.returncode}")
            if stderr:
                logger.error(f"STDERR: {stderr}")
            raise RuntimeError(f"DIAMOND blastx search failed with return code {proc.returncode}")

        if stderr:
            logger.log(99, stderr)

    except Exception as e:
        if not isinstance(e, RuntimeError):
            logger.error(f"Unexpected error in DIAMOND blastx search: {e}")
            raise RuntimeError(f"DIAMOND blastx search failed: {e}") from e
        raise

    if not Path(output_path).exists():
        logger.error(f"DIAMOND output file not found at {output_path}")
        raise RuntimeError(f"DIAMOND search did not produce output at {output_path}")
        
    logger.info(f"Successfully completed DIAMOND search. Results at {output_path}")
    return output_path, column_info_path




# ===============================
# search_filtering function
# ===============================
def search_filtering(search_output_path: str,
                     min_score_threshold: int):
    """
    Filters DIAMOND search output, keeping only the highest scoring hit per query above a minimum score threshold.

    Args:
        search_output_path (str): Path to DIAMOND search output file
        min_score_threshold (int): Minimum score threshold to keep hits

    Returns:
        search_output_path (str): Path to filtered search output file
    """
    logger.info(f"Filtering search results, keeping highest scoring hit per query (score >= {min_score_threshold})")
    
    # Determine columns structure from the file content
    columns = ["qseqid", "sseqid", "pident", "qlen", "slen", "length", "mismatch", "gapopen", "qstart", "qend",
               "sstart", "send", "evalue", "bitscore", "score"]
    
    best_hits = {}  # {qseqid: (score, full_row)}
    score_idx = columns.index('score')
    qseqid_idx = columns.index('qseqid')
    
    # Stream through file to find best hits
    with open(search_output_path, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            qseqid = fields[qseqid_idx]
            score = float(fields[score_idx])
            
            # Only keep hits with score >= threshold
            if score >= min_score_threshold:
                if qseqid not in best_hits or score > best_hits[qseqid][0]:
                    best_hits[qseqid] = (score, line)
    
    # ===========================
    # Write filtered results
    # ===========================
    temp_output = search_output_path + '.tmp'
    with open(temp_output, 'w') as f:
        for score, line in best_hits.values():
            f.write(line)
    
    # Replace original with filtered
    Path(temp_output).replace(search_output_path)
    
    logger.info(f"Filtering complete. Retained {len(best_hits)} hits after filtering.")
    return search_output_path



# =========================================
# Split search outputs per gene or dataset
# =========================================
def split_search_outputs(blast_out_path: str,
                            mapping_dict: dict,
                            output_dir: str,
                            split_column: int,
                            force: bool = False):
    """
    Splits BLAST output file into separate files per gene or dataset.

    Args:
        blast_out_path (str): Path to BLAST output file
        mapping_dict (dict): Dictionary mapping gene names or dataset names to sequence IDs
        output_dir (str): Directory to save split output files
        split_column (int): Column index in BLAST output to use for splitting
        force (bool): Whether to overwrite existing files
    
    Returns:
        split_outputs (dict): Dictionary mapping the group names to their split BLAST output file paths
    """
    split_outputs = {}

    # ====================================================================
    # Build dictionary for seqid mapping during streaming of blast output
    # ====================================================================
    seq_to_group = {}
    for group_name, seq_ids in mapping_dict.items():
        for seq_id in seq_ids:
            seq_to_group[seq_id] = group_name
    
    db_search_name = determine_dataset_name(blast_out_path, '.', 0, '').rsplit('_hits', 1)[0]
    
    logger.info(f"Processing BLAST output: {blast_out_path}")
    logger.info(f"Splitting into {len(mapping_dict)} categories") # categories = genes or datasets depending on mapping_dict
    logger.info(f"Built reverse lookup for {len(seq_to_group):,} sequence IDs")
    
    # ========================================
    # Open file handles for each split output
    # ========================================
    file_handles = {}
    for name in mapping_dict.keys():
        split_output_path = ensure_path(output_dir, f"{db_search_name}_{name}_hits.tsv", force=force)
        file_handles[name] = open(split_output_path, 'w')
        split_outputs[name] = split_output_path
    
    # =================================================================
    # Stream through BLAST output and write to appropriate split files
    # =================================================================
    line_counts = {}
    total_lines = 0
    unmatched_lines = 0
    
    try:
        with open(blast_out_path, 'r') as infile:
            for line in infile:
                total_lines += 1
                fields = line.strip().split('\t')
                
                if len(fields) > split_column:
                    split_value = fields[split_column]
                    
                    # Look up which group this belongs to
                    if split_value in seq_to_group:
                        group_name = seq_to_group[split_value]
                        file_handles[group_name].write(line)
                        line_counts[group_name] = line_counts.get(group_name, 0) + 1
                    else:
                        unmatched_lines += 1
                
                # if total_lines % 100000 == 0: 
                #     logger.debug(f"Processed {total_lines:,} lines...")

    except Exception as e:
        logger.error(f"Error while splitting BLAST output: {e}")
        raise RuntimeError(f"Failed to split BLAST output: {e}") from e
    
    # =======================
    # Close all file handles
    # =======================
    finally:
        for file_handle in file_handles.values():
            file_handle.close()
    
    logger.info(f"Total lines processed: {total_lines:,}")
    logger.info(f"Unmatched lines: {unmatched_lines:,}")
    
    for name, count in line_counts.items():
        if count > 0:
            logger.info(f"Wrote {count:,} lines for '{name}' to {split_outputs[name]}")
        else:
            logger.info(f"No hits found for '{name}'")
    
    return split_outputs





# ================================
# aastk get_hit_seqs CLI FUNCTION
# ================================
def get_hit_seqs(blast_tab: str,
                        query_path: str,
                        output_dir: str,
                        key_column: int = 0,
                        force: bool = False):
    """
    Extracts read sequences that have hits in the DIAMOND BLAST tabular output.

    Args:
        blast_tab (str): Path to DIAMOND BLAST tabular output file
        query_path (str): Path to sequencing read file, can be gzipped
        output_dir (str): Directory to save output files
        key_column (int): Column index for query sequence ID (default: 0)
        force (bool): Whether to overwrite existing files

    Returns:
        hit_seqs_path (str): Path to the extracted hit sequences in gzipped FASTQ format
        id_file (str): Path to text file containing extracted IDs
        dataset_dict (dict): Dictionary with dataset name and list of query IDs 
    """
    # Check if seqkit is available
    check_dependency_availability("seqkit")

    # ===================
    # Output file setup
    # ===================
    protein_name = determine_dataset_name(blast_tab, '.', 0, '_hits')
    out_fastq = ensure_path(output_dir, f"{protein_name}_matched.fastq.gz", force=force)
    
    # =========================================
    # Extract matching IDs from BLAST results
    # =========================================
    logger.info(f"Extracting unique sequence IDs from {blast_tab}")

    # Extract unique IDs from dereplicated BLAST results
    unique_ids = extract_unique_keys(blast_tab, key_column)

    # Fill in dataset_dict (dataset_name:[matched_seq_ids]) for downstream use
    dataset_name = determine_dataset_name(query_path, '.', 0, '')
    dataset_dict = {}
    dataset_dict[dataset_name] = list(unique_ids)

    # ======================================================
    # Write unique IDs to file for seqkit grep
    # ======================================================
    id_file = ensure_path(output_dir, f"{protein_name}_matched_ids.txt", force=force)
    with open(id_file, 'w') as f:
        for seq_id in unique_ids:
            f.write(f"{seq_id}\n")
    
    logger.info(f"Extracted {len(unique_ids)} unique sequence IDs to {id_file}")

    # ======================================================
    # Extract and write matching sequences in gzipped FASTQ
    # ======================================================
    cmd = ["seqkit", "grep", "-f", id_file, query_path, "-o", out_fastq]

    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        logger.info(f"Extracted matching sequences to {out_fastq}")
    
    except subprocess.CalledProcessError as e:
        logger.error(f"seqkit grep failed: {e.stderr}")
        raise RuntimeError(f"Failed to extract matching sequences: {e.stderr}") from e
    
    logger.info(f"Successfully extracted {len(unique_ids)} hit sequences to {out_fastq}")

    return out_fastq, id_file, dataset_dict




# ================================
# aastk merge_hits CLI FUNCTION
# ================================
def merge_hits(all_hit_seqs_paths: list,
                output_dir: str,
                threads: int,
                use_existing_merged: bool = False,
                force: bool = False):
    """
    Merges multiple FASTQ files containing hit sequences into a single FASTQ file.

    Args:
        all_hit_seqs_paths (list): List of paths to FASTQ files containing hit sequences
        output_dir (str): Directory to save merged FASTQ file
        use_existing_merged (bool): Whether to reuse existing merged file if it exists and matches the input files
        threads (int): Number of threads to use for merging
        force (bool): Whether to overwrite existing files

    Returns:
        merged_hits_path (str): Path to merged FASTQ file
        infile_list (str): Path to text file containing list of input FASTQ files used for merging
    """
    # Check if seqkit is available
    check_dependency_availability("seqkit")

    # ===================
    # Output file setup
    # ===================
    merged_hits_file = ensure_path(output_dir, "all_matched.fastq.gz", force=force)
    infile_list = ensure_path(output_dir, "all_hit_seqs_paths.txt", force=force)

    # =========================================================
    # Write input file list and check for existing merged file
    # =========================================================
    with open(infile_list, 'w') as f:
        f.write('\n'.join(all_hit_seqs_paths))

    if Path(merged_hits_file).exists() and not force:
        if use_existing_merged:
            logger.warning(
                f"Using existing merged hits file at {merged_hits_file}. "
                "Ensure it matches the current infile list."
            )
            return merged_hits_file, infile_list
        logger.error(
            f"Merged hits file already exists at {merged_hits_file}. "
            "Set use_existing_merged=True to reuse it, or force=True to rebuild."
        )
        raise FileExistsError(f"Stale merged hits file detected: {merged_hits_file}")
        
    # =======================================================
    # Merge and deduplicate all matched sequences (search 1)
    # =======================================================
    logger.info(f"Merging and deduplicating {len(all_hit_seqs_paths)} hit FASTQ files")
    
    cmd = ["seqkit", "rmdup", "-s", "-o", merged_hits_file, *all_hit_seqs_paths]

    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        logger.info(f"Successfully merged and deduplicated {len(all_hit_seqs_paths)} FASTQ files into {merged_hits_file}")
        loggger.info(f"The merged file contains {count_fastq_records(merged_hits_file):,} unique sequences")
    
    except subprocess.CalledProcessError as e:
        logger.error(f"seqkit rmdup failed: {e.stderr}")
        raise RuntimeError(f"Failed to merge and deduplicate hit sequences: {e.stderr}") from e
    
    return merged_hits_file, infile_list




# ===============================
# aastk bsr CLI FUNCTION
# ===============================
def bsr(db_search_tab: str,
        og_search_tab: str,
        column_info_path: str,
        output_dir: str,
        force: bool = False):
    """
    Computes BSR (Blast Score Ratio) using protein BLAST tab file and outgroup BLAST tab file.

    Args:
        db_search_tab (str): Path to BLAST output file from Search 1 (gene db)
        og_search_tab (str): Path to BLAST output file from Search 2 (outgroup db)
        column_info_path (str): Path to JSON file containing column index information for BLAST output
        output_dir (str): Directory to save BSR output file
        force (bool): Whether to overwrite existing files

    Returns:
        bsr_output_path (str): Path to BSR output file containing qseqid, sseqid_db, score_db, sseqid_og, score_og, bsr
    """
    # ==================
    # Output file setup
    # ==================
    protein_name = determine_dataset_name(db_search_tab, '.', 0, '_hits')
    bsr_file = ensure_path(output_dir, f"{protein_name}_bsr.tsv", force=force)

    logger.info(f"Computing blast score ratio (BSR) for {protein_name}")
    logger.info(f"Using gene blast tab file: {db_search_tab}")
    logger.info(f"Using outgroup blast tab file: {og_search_tab}")
    
    # ==========================================================
    # Load column info if provided and retrieve columns indexes
    # ==========================================================
    if column_info_path and Path(column_info_path).exists():
        try:
            with open(column_info_path, 'r') as f:
                column_info = json.load(f)
                
                # Get all required column indexes from .json file
                qseqid_col = column_info.get('qseqid', 0)
                sseqid_col = column_info.get('sseqid', 1)
                pident_col = column_info.get('pident', 2)
                length_col = column_info.get('length', 5)
                score_column = column_info.get('score', 14)
                
                logger.info(f"Loaded column indexes from {column_info_path}")
                logger.debug(f"Column indexes: qseqid={qseqid_col}, sseqid={sseqid_col}, "
                           f"pident={pident_col}, length={length_col}, score={score_column}")

        except Exception as e:
            logger.warning(f"Failed to load column info from {column_info_path}: {e}")
            # Fall back to default indexes
            qseqid_col = 0
            sseqid_col = 1
            pident_col = 2
            length_col = 5
            score_column = 14
    else:
        # Use default column indexes if no .json file provided
        qseqid_col = 0
        sseqid_col = 1
        pident_col = 2
        length_col = 5
        score_column = 14
        logger.info(f"Using default column indexes (no column_info file provided)")

    # =====================================
    # Load BLAST outputs and calculate BSR
    # =====================================
    # Read outgroup results into dict
    logger.info("Loading database search results into memory...")
    db_data = {}  # {qseqid: (sseqid_db, pident_db, length_db, score_db)}
    
    with open(db_search_tab, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            qseqid = fields[qseqid_col]
            sseqid_db = fields[sseqid_col]
            pident_db = float(fields[pident_col])
            length_db = float(fields[length_col])
            score_db = float(fields[score_column])
            db_data[qseqid] = (sseqid_db, pident_db, length_db, score_db)
    
    logger.info(f"Loaded {len(db_data):,} database hits")
    
    # Stream through outgroup results and calculate BSR
    logger.info("Streaming through outgroup search results and calculating BSR...")
    matched_count = 0
    total_count = 0
    
    with open(bsr_file, 'w') as out_f:
        # Write header
        out_f.write('qseqid\tsseqid_db\tpident_db\tlength_db\tscore_db\tsseqid_og\tpident_og\tlength_og\tscore_og\tBSR\n')
        
        with open(og_search_tab, 'r') as og_f:
            for line in og_f:
                fields = line.strip().split('\t')
                qseqid = fields[qseqid_col]
                total_count += 1
                
                # Only process if we have database data for this query
                if qseqid in db_data:
                    sseqid_og = fields[sseqid_col]
                    pident_og = float(fields[pident_col])
                    length_og = float(fields[length_col])
                    score_og = float(fields[score_column])
                    
                    sseqid_db, pident_db, length_db, score_db = db_data[qseqid]
                    
                    # Calculate BSR
                    bsr_value = score_db / score_og if score_og > 0 else 0.0
                    
                    # Write result
                    out_f.write(f"{qseqid}\t{sseqid_db}\t{pident_db}\t{length_db}\t{score_db}\t{sseqid_og}\t{pident_og}\t{length_og}\t{score_og}\t{bsr_value}\n")
                    matched_count += 1
    
    logger.info(f"BSR calculated for {matched_count:,} out of {total_count:,} outgroup hits")
    logger.info(f"BSR results written to {bsr_file}")
    
    return bsr_file




# ===============================
# aastk rasr_plot CLI FUNCTION
# ===============================
def rasr_plot(bsr_file: str,
                output_dir: str,
                bsr_cutoff: float = 0.9,
                dbmin: int = 110,
                force: bool = False):
    """
    Creates a scatter plot of the BSR data flanked by histograms showing the distribution of datapoints alongside the axes.
    
    Args:
        bsr_file (str): Path to BSR results file
        output_dir (str): Directory to save output plot
        bsr_cutoff (float): Minimum BSR threshold to draw as reference line (default: 0.9)
        dbmin (int): Minimum database score threshold to draw as reference line (default: 110)
        force (bool): Whether to overwrite existing files
    """
    logger = logging.getLogger(__name__)

    # Handle None values from cli.py
    if bsr_cutoff is None:
        bsr_cutoff = 0.9
    if dbmin is None:
        dbmin = 110

    protein_name = determine_dataset_name(bsr_file, '.', 0, '_bsr')

    # ===============================
    # Output file path setup
    # ===============================
    out_graph = ensure_path(output_dir, f'{protein_name}_bsr.png', force=force)

    logger.info(f"Creating BSR scatter plot for {protein_name}")

    try:
        # ===============================
        # Load data
        # ===============================
        bsr_df = pd.read_csv(bsr_file, sep='\t', header=0)

        required_cols = ['score_db', 'score_og', 'pident_db', 'BSR']
        for col in required_cols:
            if col not in bsr_df.columns:
                raise ValueError(f"Required column '{col}' not found in BSR data")

        # ===============================
        # Define axis ranges
        # ===============================
        x_max = bsr_df['score_og'].max()
        y_max = bsr_df['score_db'].max()
        x_min = 45
        y_min = 45

        # Add padding
        x_range = x_max - x_min
        y_range = y_max - y_min
        x_pad = max(x_range * 0.05, 10)
        y_pad = max(y_range * 0.05, 10)

        xlim = (max(0, x_min - x_pad), x_max + x_pad)
        ylim = (max(0, y_min - y_pad), y_max + y_pad)

        # ===============================
        # Layout
        # ===============================
        fig, ax = plt.subplots(figsize=(8, 8), layout='constrained')

        # ===============================
        # Scatter plot
        # ===============================
        scatter = ax.scatter(
            bsr_df['score_og'],
            bsr_df['score_db'],
            c=bsr_df['pident_db'],
            cmap='viridis',
            edgecolor='k',
            alpha=0.5,
            vmin=0,
            vmax=100
        )
        # Set labels and title
        ax.set_xlabel('Alignment score to outgroup')
        ax.set_ylabel('Alignment score to database')
        ax.set_title(f'BSR Plot for {protein_name}')
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        # Add cutoff lines
        if dbmin is not None and dbmin > 0:
            ax.axhline(dbmin, color='black', linestyle='--', linewidth=1.0, label=f'DB score ≥ {dbmin}')
        if bsr_cutoff is not None and bsr_cutoff > 0:
            x_vals = np.linspace(xlim[0], xlim[1], 500)
            y_vals = bsr_cutoff * x_vals
            ax.plot(x_vals, y_vals, color='black', linestyle='--', linewidth=1.0, label=f'BSR ≥ {bsr_cutoff}')
        
        # Add colorbar
        cb_ax = inset_axes(ax, width="5%", height="30%", loc='upper left', borderpad=1)
        cbar = fig.colorbar(scatter, cax=cb_ax)
        cbar.set_label('% seq id to database')

        # ===============================
        # Save
        # ===============================
        fig.savefig(out_graph, dpi=300)
        plt.close(fig)
        logger.info(f"BSR scatter plot saved to {out_graph}")

        return out_graph

    except Exception as e:
        logger.error(f"Error creating BSR plot: {e}")
        raise RuntimeError(f"Failed to create BSR plot: {e}") from e




# ===============================
# aastk pasr_select CLI FUNCTION
# ===============================
def rasr_select(score_cutoff: float,
                bsr_cutoff: float,
                matched_fastq: str,
                bsr_file: str,
                output_dir: str,
                force: bool = False):
    """
    Subsets RASR results based on BSR and database score thresholds.

    Args:
        score_cutoff (float): Minimum database score threshold for selection
        bsr_cutoff (float): Minimum BSR threshold for selection
        matched_fastq (str): Path to FASTQ file containing sequences that matched the gene database (Search 1)
        bsr_file (str): Path to BSR results file containing qseqid, sseqid_db, score_db, sseqid_og, score_og, bsr
        output_dir (str): Directory to save output files
        force (bool): Whether to overwrite existing files

    Returns:
        selected_fastq (str): Path to FASTQ file containing sequences that passed the selection criteria
        selected_ids_file (str): Path to text file containing IDs of selected sequences
    """
    # Check if seqkit is available
    check_dependency_availability("seqkit")

    logger.info(f"Selecting RASR hits with database score >= {score_cutoff} and BSR >= {bsr_cutoff}")

    # ===============================
    # Output file path setup
    # ===============================
    protein_name = determine_dataset_name(bsr_file, '.', 0, '_bsr')
    out_fastq = ensure_path(output_dir, f"{protein_name}_selected.fastq", force=force)
    id_file = ensure_path(output_dir, f"{protein_name}_selected_ids.txt", force=force)
    
    # ===================================
    # Extract reads IDs based on cutoffs
    # ===================================
    # Filter based on cutoffs and write selected IDs to file
    bsr_df = pd.read_csv(bsr_file, sep='\t', header=0)

    selected_ids = bsr_df[(bsr_df['score_db'] >= score_cutoff) 
                            & (bsr_df['BSR'] >= bsr_cutoff)]['qseqid'].unique()

    with open(id_file, 'w') as f:
        for seq_id in selected_ids:
            f.write(f"{seq_id}\n")

    logger.info(f"Selected {len(selected_ids)} sequences. IDs written to {id_file}")

    # ===============================
    # Extract selected sequences
    # ===============================
    try:
        cmd = ["seqkit", "grep", "-f", id_file, matched_fastq, "-o", out_fastq]

        subprocess.run(cmd, check=True)
        logger.info(f"Extracted selected sequences to {out_fastq}")

    
    except Exception as e:
        logger.error(f"seqkit grep failed: {e.stderr}")
        raise RuntimeError(f"Failed to extract selected sequences: {e.stderr}") from e

    return out_fastq, id_file




# ===============================
# aastk rasr WORKFLOW
# ===============================
def rasr_multiple(query: str,
            gene_db: str,
            outgrp_db: str,
            output_dir: str,
            sensitivity: str,
            block: int,
            chunk: int,
            threads: int = 1,
            force: bool = False,
            keep: bool = False,
            key_column: int = 0,
            bsr_cutoff: float = 0.9,
            dbmin: int = 110,
            use_existing_merged: bool = False):
    """
    RASR workflow:
    Runs:
        Once:
            aastk list_inputs       # Create lists of input files for datasets and gene dbs
            aastk merge_dbs         # Merge gene dbs in one uniq FASTA if multiple
            aastk build             # Build diamond db of genes if not exists

        For each dataset:
            aastk search            # Search merged database (Search 1) (DS1_hits.tsv)
            aastk split             # Seperate Search 1 results per gene (DS1_GENE1.tsv, DS1_GENE2.tsv, etc)
            aastk get_hit_seqs      # Extract hits for each gene (DS1_GENE1_matched.fastq, DS1_GENE2_matched.fastq, etc)
        
        Once:
            aastk merge_hits        # Merge all extracted hits from all datasets, all genes (all_hits.fastq)
            aastk search            # Search outgroup db with merged hits (Search 2) (og_allhits.tsv)
            aastk split             # Seperate Search 2 results per dataset (og_DS1.tsv, og_DS2.tsv, etc)
        
        For each dataset: 
            For each gene:
                aastk bsr               # Calculate BSR values from Search 1 and Search 2 for each gene (DS1_GENE1_bsr.tsv, DS1_GENE2_bsr.tsv, etc)
                aastk rasr_plot         # Plot BSR values for each gene and dataset (DS1_GENE1_bsr_plot.png, DS1_GENE2_bsr_plot.png, etc)
                aastk select            # Select hits above BSR and dbmin thresholds (DS1_GENE1_selected.tsv, DS1_GENE2_selected.tsv, etc)
                

    Args:
        query (str): path to sequencing read file, can be gzipped
        gene_db (str): path to gene of interest diamond database file
        outgrp_db (str): path to outgroup diamond database file
        output_dir (str): directory to save output files
        sensitivity (str): sensitivity setting for diamond search
        block (int): block size parameter for diamond search
        chunk (int): chunk size parameter for diamond search
        threads (int): number of threads to use
        force (bool): whether to force overwrite existing files
        keep (bool): if True, keep intermediate files
        key_column (int): column index in BLAST output for query ID
        bsr_cutoff (float): minimum BSR threshold (default: 0.9)
        dbmin (int): minimum database score threshold (default: 110)
        use_existing_merged (bool): whether to reuse existing merged hits file if it exists and matches the input files
    
    Returns:
        results (dict): dictionary containing paths to all output files
    """
    # ==========================
    # Normalize inputs to lists
    # ==========================
    logger.info("Normalizing inputs to lists")

    db_files, query_files = list_inputs(gene_db, query)
    
    # ========================
    # Output directory setup
    # ========================

    logger.info(f"Running RASR workflow for {len(query_files)} dataset(s) and {len(db_files)} gene database(s)")
    logger.info(f"Output directory: {output_dir}")

    # store all the output paths in dictionaries
    intermediate_results = {}
    results_by_ds_gene = {} # {(dataset_name, gene_name): {result_info_dict}}

    try:
        # ==================================
        # Gene of interest database merging
        # ==================================
        logger.info("Merging protein databases")
        merged_gene_db_path, db_seqid_dict = merge_dbs(db_files, output_dir, force=force)

        intermediate_results['merged_gene_db_fasta'] = merged_gene_db_path

        # ===================================
        # Gene of interest database building
        # ===================================
        logger.info("Building diamond merged protein database")
        db_path = pasr_build(merged_gene_db_path, threads, output_dir, force=force)

        intermediate_results['merged_gene_db_diamond'] = db_path

        # ===================================
        # Track outputs per dataset and gene
        # ===================================
        all_hit_seqs_paths = []  # Accumulate all _matched.fastq.gz file paths across datasets and genes
        dataset_to_queries = {}  # {dataset_name: [qseqids]}

        # ==================================================
        # Step 1: Search each dataset against gene database
        # ==================================================     
        logger.info("=== STEP 1: Searching gene database ===")
        dataset_output_dirs = set()
        for query_file in query_files:
            dataset_name = determine_dataset_name(query_file, '.', 0, '')
            dataset_output_dir = ensure_path(output_dir, dataset_name, force=force)
            dataset_output_dirs.add(dataset_output_dir)
            
            logger.info(f"Processing dataset: {dataset_name}")

            # ======================
            # Search gene database
            # ======================
            logger.info(f"Searching protein database")
            search_output, column_info_path = search(db_path, str(query_file), threads, dataset_output_dir, sensitivity, block=block, chunk=chunk, force=force)
            
            # =========================
            # Filter search results
            # =========================
            search_output = search_filtering(search_output, min_score_threshold=50)

            if 'search_output' not in intermediate_results:
                intermediate_results['search_output'] = []
            intermediate_results['search_output'].append(search_output)

            # ==============================
            # Split search results per gene
            # ==============================
            logger.info(f"Splitting search results per gene for {dataset_name}")
            
            db_split_outputs = split_search_outputs(
                search_output,
                mapping_dict=db_seqid_dict,
                output_dir=dataset_output_dir,
                split_column=1,
                force=force
            )

            for gene_name, split_output in list(db_split_outputs.items()):
                gene_output_dir = ensure_path(dataset_output_dir, gene_name, force=force)
                dest_path = ensure_path(gene_output_dir, Path(split_output).name, force=force)
                Path(split_output).replace(dest_path)
                db_split_outputs[gene_name] = dest_path

            if 'db_split_outputs' not in intermediate_results:
                intermediate_results['db_split_outputs'] = []
            intermediate_results['db_split_outputs'].append(db_split_outputs)            

            # ===============================
            # Extract hit sequences per gene
            # ===============================
            logger.info(f"Extracting hit sequences per gene for {dataset_name}")
            
            for gene_name, split_output in db_split_outputs.items():
                gene_output_dir = ensure_path(dataset_output_dir, gene_name, force=force) 

                hit_seqs_path, id_file, dataset_dict = get_hit_seqs(split_output, str(query_file), 
                                                                        output_dir=gene_output_dir, 
                                                                        key_column=key_column, 
                                                                        force=force)

                all_hit_seqs_paths.append(str(hit_seqs_path))

                results_by_ds_gene[(dataset_name, gene_name)] = {
                    'dataset_name': dataset_name,
                    'gene_name': gene_name,
                    'db_search_output': split_output,
                    'matched_fastq_path': hit_seqs_path,
                    'matched_ids_path': id_file,
                    'output_dir': gene_output_dir
                }
                
                # Accumulate dataset matched queries for later outgroup splitting
                if dataset_name not in dataset_to_queries:
                    dataset_to_queries[dataset_name] = []
                dataset_to_queries[dataset_name].extend(dataset_dict[dataset_name])

        if len(dataset_output_dirs) != len(query_files):
            raise RuntimeError(
                "Dataset output directory count does not match input datasets: "
                f"{len(dataset_output_dirs)} dirs for {len(query_files)} datasets"
            )

        intermediate_results['column_info_path'] = column_info_path

        # ============================================================================
        # Step 2: Merge hits from all datasets and genes and search outgroup database
        # ============================================================================
        logger.info("=== STEP 2: Searching outgroup database ===")

        logger.info("Merging hits from all datasets and genes")
        merged_hits_file, infile_list = merge_hits(all_hit_seqs_paths, output_dir, threads, use_existing_merged=use_existing_merged, force=force)

        intermediate_results['merged_hits_file'] = merged_hits_file
        intermediate_results['infile_list'] = infile_list

        # =================== Search outgroup database with merged hits ===================

        logger.info("Searching outgroup database with merged hits")
        og_search_output, og_column_info_path = search(outgrp_db, merged_hits_file, threads, output_dir, sensitivity, block=block, chunk=chunk, force=force)

        intermediate_results['og_search_output'] = og_search_output
        intermediate_results['og_column_info_path'] = og_column_info_path

        # =================== Split outgroup search results per dataset ===================

        logger.info("Splitting outgroup search results per dataset")
        og_split_outputs = split_search_outputs(og_search_output, mapping_dict=dataset_to_queries, output_dir=output_dir, split_column=0, force=force)

        intermediate_results['og_split_outputs'] = og_split_outputs

        # ================================================================
        # Step 3: Calculate BSR and generate plots for each dataset x gene
        # ================================================================
        logger.info("=== STEP 3: Calculating BSR and generating visualizations ===")

        for (dataset_name, gene_name), result_info in results_by_ds_gene.items():
            
            gene_output_dir = result_info['output_dir']
            db_split_output = result_info['db_search_output']
            og_split_output = og_split_outputs.get(dataset_name)
            
            if db_split_output and og_split_output:
                logger.info(f"Processing {dataset_name} x {gene_name}")

                # =================== Calculate BSR values ===================
                logger.info(f"Calculate BSR values")
                bsr_output = bsr(db_split_output, og_split_output, column_info_path, output_dir=gene_output_dir, force=force)

                result_info['bsr_output'] = bsr_output

                # =================== Plot BSR values ===================

                logger.info(f"Plotting BSR values")
                plot_output = rasr_plot(bsr_output, output_dir=gene_output_dir, bsr_cutoff=bsr_cutoff, dbmin=dbmin, force=force)

                result_info['bsr_plot'] = plot_output
          
                # =================== Select hits above BSR and dbmin thresholds ===================

                logger.info(f"Selecting hits above thresholds")
                hit_seqs_path = Path(gene_output_dir) / f"{determine_dataset_name(db_split_output, '.', 0, '_hits')}_matched.fastq.gz"
                selected_fastq, selected_ids_file = rasr_select(dbmin, bsr_cutoff, str(hit_seqs_path), bsr_output, output_dir=gene_output_dir, force=force)

                result_info['selected_fastq'] = selected_fastq
                result_info['selected_ids_file'] = selected_ids_file
            else:
                logger.warning(f"Missing split outputs for {dataset_name} x {gene_name}, skipping BSR calculation")
        
        logger.info("RASR workflow completed successfully")
        return results_by_ds_gene

    except Exception as e:
        logger.error(f"Error during RASR workflow: {e}")
        raise RuntimeError(f"RASR workflow failed: {e}") from e            

    # ===============================
    # Cleanup intermediate files if keep=False
    # ===============================
    finally:
        if not keep and intermediate_results:
            logger.info("Cleaning up intermediate files")
            for key, filepath in intermediate_results.items():
                try:
                    file_path = Path(filepath)
                    if file_path.exists():
                        file_path.unlink()
                        logger.debug(f"Deleted {filepath}")
                except Exception as e:
                    logger.warning(f"Failed to delete {filepath}: {e}")
