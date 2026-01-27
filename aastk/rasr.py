#!/usr/bin/env python3
import json
from .util import *
from .pasr import build as pasr_build
import subprocess
import logging
from pathlib import Path


logger = logging.getLogger(__name__)




# ===============================
# aastk build CLI FUNCTION
# ===============================
def rasr_build(gene_db_fasta: str,
            gene_db_dir: str,
            threads: int,
            force: bool = False):
    """
    Builds a DIAMOND protein database from a fasta file.
    
    Args:
        gene_db_fasta (str): Path to gene of interest fasta file
        gene_db_dir (str): Directory where the database will be created
        threads (int): Number of threads to use
        force (bool): If True, overwrite existing database
    
    Returns:
        gene_db_path (str): Path to the created DIAMOND database

    Raises:
        RuntimeError: If the DIAMOND database creation fails.
    """
    
    # Check if diamond is available
    check_dependency_availability("diamond")

    # Check if gene_db_fasta exists 
    if Path(gene_db_fasta).is_file():
        pass
    else:
        logger.error("Input seed FASTA not found")
        raise FileNotFoundError(f"Seed FASTA file does not exist: {gene_db_fasta}")
    
    gene_db_fasta_filename = Path(gene_db_fasta).stem
    
    # Create output directory if it doesn't exist
    gene_db_out_path = ensure_path(gene_db_dir, f"{gene_db_fasta_filename}_db", force=force)

    # Log the path
    logger.info(f"Building DIAMOND database for {gene_db_fasta_filename} at {gene_db_out_path}")
 
    # =======================================
    # DIAMOND makedb command construction
    # =======================================
    cmd = ["diamond", "makedb",
           "--in", gene_db_fasta,
           "-d", gene_db_out_path,
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

        if proc.returncode != 0:
            logger.error(f"DIAMOND makedb failed with return code {proc.returncode}")
            if stderr:
                logger.error(f"STDERR: {stderr}")
            raise RuntimeError(f"DIAMOND database creation failed with return code {proc.returncode}")

        if stderr:
            logger.log(99, stderr)

    except Exception as e:
        if not isinstance(e, RuntimeError):
            logger.error(f"Unexpected error in building the DIAMOND database: {e}")
            raise RuntimeError(f"DIAMOND database creation failed: {e}") from e
        raise

    db_file = Path(f"{gene_db_out_path}.dmnd")
    if not db_file.exists():
        logger.error(f"DIAMOND database file not found at {db_file}")
        raise RuntimeError(f"DIAMOND database was not created at {db_file}")


    logger.info(f"Successfully built DIAMOND database at {gene_db_out_path}")
    return gene_db_out_path




# ===============================
# aastk search CLI FUNCTION
# ===============================
def rasr_search(gene_db_out_path: str,
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
    gene_db_filename = determine_dataset_name(gene_db_out_path, '.', 0, '_db')

    # ===============================
    # Output file path setup
    # ===============================
    output_path = ensure_path(output_dir, f"{gene_db_filename}_hits.txt", force=force)
    column_info_path = ensure_path(output_dir, f"{gene_db_filename}_columns.json", force=force)

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
    cmd = ["diamond", "blastx",
           "-d", gene_db_out_path,
           "-q", query_fastq,
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
def rasr_get_hit_seqs(blast_tab: str,
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
# aastk rasr WORKFLOW
# ===============================
def rasr(query: str,
            gene_db_fasta: str,
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
        query_fastq (str): path to sequencing read file, can be gzipped
        gene_db_fasta (str): path to gene of interest diamond database file
        output_dir (str): directory to save output files
        sensitivity (str): sensitivity setting for diamond search
        block (int): block size parameter for diamond search
        chunk (int): chunk size parameter for diamond search
        threads (int): number of threads to use
        force (bool): whether to force overwrite existing files
        keep (bool): if True, keep intermediate files; if False, delete them after workflow completion
        key_column (int): column index in BLAST output to use as key for sequence extraction
    
    """

    protein_name = Path(gene_db_fasta).stem

    # check for output_dir
    logger.info(f"Running RASR workflow for {protein_name}")
    logger.info(f"Output directory: {output_dir}")

    # store all the output paths in dictionaries
    intermediate_results = {}
    results = {}

    try:
        # ===============================
        # Gene of interest database building
        # ===============================
        logger.info("Building protein database")
        db_path = pasr_build(gene_db_fasta, threads, output_dir, force=force)

        # intermediate_results['db_path'] = f"{db_path}.dmnd"
        results['db_path'] = f"{db_path}.dmnd"

        # ===============================
        # Gene of interest database search
        # ===============================
        logger.info("Searching protein database")

        search_output, column_info_path = rasr_search(db_path, query, threads, output_dir, sensitivity, block, chunk, force=force)

        results['search_output'] = search_output
        results['column_info_path'] = column_info_path

        # ===============================
        # Read sequence extraction
        # ===============================
        logger.info("Extracting hit sequences from search results")
        matched_fastq, matched_stats = rasr_get_hit_seqs(search_output, query, output_dir, key_column=0, force=force)

        results['hit_seqs_path'] = matched_fastq
        results['hit_seqs_stats'] = matched_stats

        # ===============================
        # Outgroup database search
        # ===============================

        # ===============================
        # Score calculations
        # ===============================

        # ===============================
        # Visualization
        # ===============================

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