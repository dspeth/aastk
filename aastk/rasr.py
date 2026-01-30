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
    output_path = ensure_path(output_dir, f"{db_filename}_hits.txt", force=force)
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
           "--gapopen", str(14), 
           "--gapextend", str(2),
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
# aastk get_hit_seqs CLI FUNCTION
# ===============================
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
        force (bool): Whether to overwrite existing files

    Returns:
        hit_seqs_path (str): Path to the extracted hit sequences in FASTA format
    """
    # Check if seqkit is available
    check_dependency_availability("seqkit")

    # ===================
    # Output file setup
    # ===================
    protein_name = determine_dataset_name(blast_tab, '.', 0, '_hits')

    id_file = ensure_path(output_dir, f"{protein_name}_ids.txt", force=force)
    out_fastq = ensure_path(output_dir, f"{protein_name}_matched.fastq", force=force)
    stats_path = ensure_path(output_dir, f"{protein_name}_matched.stats", force=force)

    # =============================================
    # Extract matching reads IDs from BLAST results
    # =============================================
    extract_cmd = f"cut -f{key_column + 1} {blast_tab} | sort -u > {id_file}"
    subprocess.run(extract_cmd, shell=True, check=True)

    logger.info(f"Extracted unique sequence IDs to {id_file}")

    # =====================================
    # Extract and write matching sequences
    # =====================================
    cmd = ["seqkit", "grep", "-f", id_file, query_path, "-o", out_fastq]

    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        logger.info(f"Extracted matching sequences to {out_fastq}")
    
    except subprocess.CalledProcessError as e:
        logger.error(f"seqkit grep failed: {e.stderr}")
        raise RuntimeError(f"Failed to extract matching sequences: {e.stderr}") from e
    
    logger.info(f"Successfully extracted hit sequences to {out_fastq}")

    # ===============================
    # Generate sequence statistics
    # ===============================
    cmd = ["seqkit", "stats",
            out_fastq,
            "-o", stats_path,
            "-a"]

    logger.debug(f"Running command: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

    return out_fastq, stats_path




# ===============================
# aastk bsr CLI FUNCTION
# ===============================
def bsr(gene_blast_tab: str,
            outgrp_blast_tab: str,
            output_dir: str,
            column_info_path: str = None,
            force: bool = False):
    """
    Computes BSR (Blast Score Ratio) using protein BLAST tab file and outgroup BLAST tab file.
    Args:
        gene_blast_tab (str): Path to protein BLAST tab file against gene database
        outgrp_blast_tab (str): Path to protein BLAST tab file against outgroup database
        output_dir (str): Directory to save output files
        column_info_path (str): Path to JSON file with column information
        force (bool): Whether to overwrite existing files
    Returns:
        bsr_file (str): Path to the output file with BSR values.
    """
    # Validate input files
    if not Path(gene_blast_tab).exists():
        logger.error(f"Gene BLAST file not found: {gene_blast_tab}")
        raise FileNotFoundError(f"Gene BLAST file does not exist: {gene_blast_tab}")
    if not Path(outgrp_blast_tab).exists():
        logger.error(f"Outgroup BLAST file not found: {outgrp_blast_tab}")
        raise FileNotFoundError(f"Outgroup BLAST file does not exist: {outgrp_blast_tab}")

    # ===============================
    # Output file setup
    # ===============================
    protein_name = determine_dataset_name(gene_blast_tab, '.', 0, '_hits')

    bsr_file = ensure_path(output_dir, f"{protein_name}_bsr.tsv", force=force)

    logger.info(f"Computing blast score ratio (BSR) for {protein_name}")
    logger.info(f"Using gene blast tab file: {gene_blast_tab}")
    logger.info(f"Using outgroup blast tab file: {outgrp_blast_tab}")

    # Load column info if provided and retrieve columns indices
    if column_info_path and Path(column_info_path).exists():
        try:
            with open(column_info_path, 'r') as f:
                column_info = json.load(f)
                
                # Get all required column indices from JSON
                qseqid_col = column_info.get('qseqid', 0)
                sseqid_col = column_info.get('sseqid', 1)
                pident_col = column_info.get('pident', 2)
                score_column = column_info.get('score', 14)
                
                logger.info(f"Loaded column indices from {column_info_path}")
                logger.debug(f"Column indices: qseqid={qseqid_col}, sseqid={sseqid_col}, "
                           f"pident={pident_col}, score={score_column}")

        except Exception as e:
            logger.warning(f"Failed to load column info from {column_info_path}: {e}")
            # Fall back to default indices
            qseqid_col = 0
            sseqid_col = 1
            pident_col = 2
            score_column = 14
    else:
        # Use default column indices if no JSON file provided
        qseqid_col = 0
        sseqid_col = 1
        pident_col = 2
        score_column = 14
        logger.info(f"Using default column indices (no column_info file provided)")
        
    
    # Load blast results into DataFrames
    dtype={0: "string", 1: "string", 2: "float32", 3: "float32",}
    usecols = [qseqid_col, sseqid_col, pident_col, score_column]

    gene_df = pd.read_csv(gene_blast_tab, sep='\t', usecols=usecols, dtype=dtype,
                            names=['qseqid', 'sseqid_db', 'pident_db', 'score_db'])
    outgrp_df = pd.read_csv(outgrp_blast_tab, sep='\t', usecols=usecols, dtype=dtype,
                           names=['qseqid', 'sseqid_og', 'pident_og', 'score_og'])

    # Merge DataFrames on 'qseqid'
    bsr_df = gene_df.merge(outgrp_df, on='qseqid', how='inner') # inner join to keep only matching qseqid (default behavior)

    # Calculate BSR
    bsr_df["BSR"] = np.where(bsr_df["score_og"] > 0, bsr_df["score_db"] / bsr_df["score_og"], 0.0)

    # Save BSR results to file
    bsr_df.to_csv(bsr_file, sep='\t', index=False)
    logger.info(f"BSR results written to {bsr_file}")

    return bsr_file




# ===============================
# aastk rasr_plot CLI FUNCTION
# ===============================
def rasr_plot(bsr_file: str,
             output_dir: str,
             force: bool = False):
    """
    Creates a scatter plot of the BSR data.
    """
    logger = logging.getLogger(__name__)

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

        required_cols = ['score_db', 'score_og', 'pident_db']
        for col in required_cols:
            if col not in bsr_df.columns:
                raise ValueError(f"Required column '{col}' not found in BSR data")
        
        fig, ax = plt.subplots(figsize=(8, 8))

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
        ax.set_xlabel('Alignment score to outgroup')
        ax.set_ylabel('Alignment score to database')
        ax.set_title(f'BSR Plot for {protein_name}')
        
        # Add colorbar
        cbar = fig.colorbar(scatter, ax=ax)
        cbar.set_label('% seq id to database')
        
        # Adjust layout
        plt.tight_layout()

        # ===============================
        # Save
        # ===============================
        fig.savefig(out_graph, dpi=300)
        plt.close(fig)

        logger.info(f"Successfully created BSR plot at {out_graph}")
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
    """
    # Check if seqkit is available
    check_dependency_availability("seqkit")

    logger.info(f"Selecting RASR hits with score >= {score_cutoff} and BSR >= {bsr_cutoff}")

    # ===============================
    # Output file path setup
    # ===============================
    protein_name = determine_dataset_name(bsr_file, '.', 0, '_bsr')
    out_fastq = ensure_path(output_dir, f"{protein_name}_selected.fastq", force=force)
    id_file = ensure_path(output_dir, f"{protein_name}_selected_ids.txt", force=force)
    stats_file = ensure_path(output_dir, f"{protein_name}_selected.stats", force=force)

    
    # ===============================
    # Load BSR data
    # ===============================
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

        # ===============================
        # Generate statistics
        # ===============================
        cmd = ["seqkit", "stats",
                out_fastq,
                "-o", stats_file,
                "-a"]

        subprocess.run(cmd, check=True)
        logger.info(f"Generated statistics for selected sequences at {stats_file}")
        return out_fastq, stats_file
    
    except Exception as e:
        logger.error(f"seqkit grep failed: {e.stderr}")
        raise RuntimeError(f"Failed to extract selected sequences: {e.stderr}") from e

        


# ===============================
# aastk rasr WORKFLOW
# ===============================
def rasr(query: str,
            gene_db_fasta: str,
            outgrp_db: str,
            output_dir: str,
            sensitivity: str,
            block: int,
            chunk: int,
            threads: int = 1,
            force: bool = False,
            keep: bool = False,
            key_column: int = 0):
    """
    RASR workflow:
    Runs:
        aastk build
        aastk search    
        aastk get_hit_seqs  # catch reads that hit gene db
        aastk search        # against_outgroup
        aastk bsr
        aastk rasr_plot

    Args:
        query (str): path to sequencing read file, can be gzipped
        gene_db_fasta (str): path to gene of interest diamond database file
        outgrp_db (str): path to outgroup diamond database file
        output_dir (str): directory to save output files
        sensitivity (str): sensitivity setting for diamond search
        block (int): block size parameter for diamond search
        chunk (int): chunk size parameter for diamond search
        threads (int): number of threads to use
        force (bool): whether to force overwrite existing files
        keep (bool): if True, keep intermediate files; if False, delete them after workflow completion
        key_column (int): column index in BLAST output to use as key for sequence extraction
    
    Returns:
        results (dict): dictionary containing paths to all output files
    """
    # ==============================
    # Output directory setup
    # ==============================
    protein_name = Path(gene_db_fasta).stem
    dataset_name = determine_dataset_name(query, '.', 0, '')

    protein_dir = ensure_path(output_dir, f"rasr_out_{protein_name}", force=force)
    dataset_output_dir = ensure_path(protein_dir, dataset_name, force=force)

    logger.info(f"Running RASR workflow for {protein_name}")
    logger.info(f"Query dataset: {dataset_name}")
    logger.info(f"Output directory: {output_dir}")

    # store all the output paths in dictionaries
    intermediate_results = {}
    results = {}

    try:
        # ===============================
        # Gene of interest database building
        # ===============================
        logger.info("Building protein database")
        db_path = pasr_build(gene_db_fasta, threads, protein_dir, force=force)

        # intermediate_results['db_path'] = f"{db_path}.dmnd"
        results['db_path'] = f"{db_path}.dmnd"

        # ===============================
        # Gene of interest database search
        # ===============================
        logger.info("Searching protein database")

        search_output, column_info_path = search(db_path, query, threads, dataset_output_dir, sensitivity, block, chunk, force=force)

        results['search_output'] = search_output
        results['column_info_path'] = column_info_path

        # ===============================
        # Read sequence extraction
        # ===============================
        logger.info("Extracting hit sequences from search results")
        matched_fastq, matched_stats = get_hit_seqs(search_output, query, dataset_output_dir, key_column=key_column, force=force)
        results['hit_seqs_path'] = matched_fastq
        results['hit_seqs_stats'] = matched_stats

        # ===============================
        # Outgroup database search
        # ===============================
        logger.info("Searching outgroup database")
        outgrp_search_output, outgrp_column_info_path = search(outgrp_db, matched_fastq, threads, dataset_output_dir, sensitivity, block, chunk, force=force)

        results['outgrp_search_output'] = outgrp_search_output
        results['outgrp_column_info_path'] = outgrp_column_info_path

        # ===============================
        # Score calculations
        # ===============================
        logger.info("Calculating BSR values")
        bsr_output = bsr(search_output, outgrp_search_output, dataset_output_dir, column_info_path=column_info_path, force=force)
        
        results['bsr_output'] = bsr_output

        # ===============================
        # Visualization
        # ===============================
        logger.info("Creating BSR plot")
        plot_output = rasr_plot(bsr_output, dataset_output_dir, force=force)

        results['plot_output'] = plot_output

        # ===============================
        # Final selection of RASR hits
        # ===============================
        logger.info("Selecting final RASR hits")
        selected_fastq, selected_stats = rasr_select(score_cutoff=100.0, bsr_cutoff=0.9, matched_fastq=matched_fastq, bsr_file=bsr_output, output_dir=dataset_output_dir, force=force)

        logger.info("RASR workflow completed successfully")
        return results

    except Exception as e:
        logger.error(f"RASR workflow failed: {e}")
        exit()
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