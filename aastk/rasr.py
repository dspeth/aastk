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
# aastk build_execution_plan FUNCTION
# ===============================
def build_execution_plan(query,
                         query_dir,
                         seed,
                         seed_dir):
    """
    Builds an execution plan for the RASR pipeline based on the provided input combinations of query and seed files.

    Args:
        query (str): Path to a single FASTQ file containing sequencing reads to search against the seed database.
        query_dir (str): Path to a directory containing multiple FASTQ files to search against the seed database.
        seed (str): Path to a single FASTA file containing reference sequences for the gene of interest to build the seed database.
        seed_dir (str): Path to a directory containing multiple FASTA files, each representing reference sequences for a different gene of interest to build the seed database.
    
    Returns:
        execution_plan (dict): A dictionary containing the execution plan for the RASR pipeline,
            including which steps to run based on the input combinations and the relevant information for each input type.
    """
    # ===============================================================
    # Build query_info and seed_info dictionaries for downstream use
    # ===============================================================
    # Build query_info{"path": Path, "mode": "single" or "batch", "label": "query", "expected_type": "fastq"}
    if (query is None) == (query_dir is None):
        raise ValueError("Invalid input combination. Use exactly one of --query or --query_dir.")
    
    if query is not None:
        query_path = Path(query)
        if not query_path.exists():
            raise FileNotFoundError(f"Query file does not exist: {query}")
        if not query_path.is_file():
            raise ValueError(f"Query path is not a file: {query}. --query expects a path to a single FASTQ file.")

        query_info = {"path": query_path, "mode": "single", "label": "query", "expected_type": "fastq"}

    else:
        query_path = Path(query_dir)
        if not query_path.exists():
            raise FileNotFoundError(f"Query directory does not exist: {query_dir}")
        if not query_path.is_dir():
            raise ValueError(f"Query path is not a directory: {query_dir}. --query_dir expects a path to a directory containing FASTQ files.")

        query_info = {"path": query_path, "mode": "batch", "label": "query", "expected_type": "fastq"}

    # Build seed_info{"path": Path, "mode": "single" or "batch", "label": "seed", "expected_type": "fasta"}
    if (seed is None) == (seed_dir is None):
        raise ValueError("Invalid input combination. Use exactly one of --seed or --seed_dir.")
    
    if seed is not None:
        seed_path = Path(seed)
        if not seed_path.exists():
            raise FileNotFoundError(f"Seed file does not exist: {seed}")
        if not seed_path.is_file():
            raise ValueError(f"Seed path is not a file: {seed}. --seed expects a path to a single FASTA file.")

        seed_info = {"path": seed_path, "mode": "single", "label": "seed", "expected_type": "fasta"}

    else:
        seed_path = Path(seed_dir)
        if not seed_path.exists():
            raise FileNotFoundError(f"Seed directory does not exist: {seed_dir}")
        if not seed_path.is_dir():
            raise ValueError(f"Seed path is not a directory: {seed_dir}. --seed_dir expects a path to a directory containing FASTA files.")

        seed_info = {"path": seed_path, "mode": "batch", "label": "seed", "expected_type": "fasta"}

    # ===============================
    # Build execution plan dict of steps to run based on input combinations
    # ===============================
    # The execution plan is a dictionary that specifies which steps of the pipeline to run based on the input combinations.

    seed_is_multi = seed_info["mode"] == "batch"
    query_is_multi = query_info["mode"] == "batch"

    execution_plan = {
        "query_info": query_info,
        "seed_info": seed_info,
        "steps": {
            "merge_seed_databases": seed_is_multi,                # Only merge if multiple seed files
            "build_db": True,                          # Always build database from seed files
            "search_1": True,                          # Always run search
            "filter_search": True,                     # Always run filtering
            "parse_search_1_per_db": seed_is_multi,    # Only split if multiple seed files
            "extract_hit_seqs": True,                      # Always get hit sequences
            "merge_extracted_hit_seqs": (seed_is_multi or query_is_multi),  # Merge when multiple hit files are expected
            "search_2": True,                          # Always run outgroup search
            "parse_search_2_per_query": query_is_multi,   # Only split if multiple query files
            "bsr": True,                                  # Always calculate BSR
            "rasr_bsr_plot": True,                            # Always create plot
            "rasr_select": True                           # Always select hits based on BSR and db score cutoffs
        }
    }

    # Log execution plan for debugging
    logger.info(f"Query input: {query_info['path']} (mode: {query_info['mode']})")
    logger.info(f"Seed input: {seed_info['path']} (mode: {seed_info['mode']})")

    for step, state in execution_plan["steps"].items():
        logger.info(f"  {step}: {'RUN' if state else 'SKIP'}")

    return execution_plan




# ===============================
# aastk list_inputs FUNCTION
# ===============================
def list_inputs(input_info: dict):
    """
    Lists input files path(s) based on the provided input information dictionary.

    Args:
        input_info (dict): Dictionary containing "path", "mode" ("single" | "batch"), "label"("query" | "seed"), and "expected_type"("fastq" | "fasta") for the input

    Returns:
        list_of_paths (list): List of Path objects for input files
    """
    # =========================
    # Set up input's variables
    # =========================
    input_path = input_info["path"]
    expected_type = input_info["expected_type"]
    input_mode = input_info["mode"]
    label = input_info["label"]

    list_of_paths = []

    # =================================
    # List input file(s) based on mode
    # =================================
    if input_mode == "batch":

        for element in input_path.iterdir():
            if not element.is_file():
                continue

            file_type = determine_file_type(element)

            if file_type == expected_type:
                list_of_paths.append(element)
            else:
                logger.warning(f"Skipping {element} because it is not a {expected_type.upper()} file.")
                continue

        list_of_paths.sort()
    
    elif input_mode == "single":

        file_type = determine_file_type(input_path)

        if file_type != expected_type:
            raise ValueError(f"{label} input file has invalid type: {input_path}. Expected a {expected_type.upper()} file.")

        list_of_paths.append(input_path)

    logger.info(
        f"[{label.upper()} INPUT] label={label} mode={input_mode} expected_type={expected_type} "
        f"files={len(list_of_paths)} path={input_path}"
    )

    if label == "query" and input_mode == "batch" and len(list_of_paths) > 20:
        logger.warning(
            f"Found {len(list_of_paths)} query files in batch mode at {str(input_path)}. "
            "This may lead to long runtimes and high resource usage. "
            "Consider using single mode or reducing the number of files per batch."
        )

    return list_of_paths




# ==============================
# merge_seed_databases function
# ==============================
def merge_seed_databases(db_files: list,
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
    merged_fasta_path = ensure_path(output_dir, "merged_seed_db.fasta", force=force)
    
    if Path(merged_fasta_path).exists() and not force:
        logger.info(f"Merged database already exists at {merged_fasta_path}, using existing.")
        
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
    
    logger.info(f"[MERGE_DB_DONE] Merged database created at {merged_fasta_path}")
    
    return str(merged_fasta_path), db_dict



# ================
# search function
# ================
def search(seed_db_out_path: str,
            query_fastq: str,
            threads: int,
            output_dir: str,
            sensitivity: str,
            bit_score_cutoff: int,
            max_target_seqs: int,
            block: int = 6,
            chunk: int = 2,
            force: bool = False):
    """
    Searches a DIAMOND reference database for homologous sequences.
    Single database vs. single query

    Args:
        seed_db_out_path (str): Path to DIAMOND protein database for the gene of interest
        query_fastq (str): Path to sequencing read file, can be gzipped
        threads (int): Number of threads to use
        output_dir (str): Directory to save output files
        sensitivity (str): Sensitivity setting for DIAMOND search
        bit_score_cutoff (int): Minimum bit score threshold for filtering search results
        block (int): Block size parameter for DIAMOND search
        chunk (int): Chunk size parameter for DIAMOND search
        force (bool): Whether to overwrite existing files

    Returns:
        blast_out_path: Path to tabular BLAST output file.
    """
    
    # Check if diamond is available
    check_dependency_availability("diamond")

    # automatically name files
    db_filename = determine_dataset_name(seed_db_out_path, '.', 0, '_db')

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
               "sstart", "send", "evalue", "bitscore", "score", "qtitle"]

    # save column information to json file as to not hardcode the score column in the bsr function
    column_info = {col: idx for idx, col in enumerate(columns)}
    with open(column_info_path, 'w') as f:
        json.dump(column_info, f, indent=2)

    # ==================
    # Parameter logging
    # ==================
    logger.info(
        f"[SEARCH] db={seed_db_out_path} query={query_fastq} out={output_path} col_info={column_info_path} \n",
        f"sensitivity={sensitivity if sensitivity else 'default'} block={block} chunk={chunk} min_score={bit_score_cutoff} max_target_seqs={max_target_seqs} algo=0"
    )

    # ====================================
    # DIAMOND blastx command construction
    # ====================================
    cmd = ["diamond", "blastx",
           "-d", seed_db_out_path,
           "-q", query_fastq,
           "-p", str(threads),
           "-o", output_path,
           "--max-target-seqs", str(max_target_seqs),
           "--matrix", "blosum45",
           "--masking", str(0),
           "--outfmt", str(6), *columns,
           "-b", str(block),
           "-c", str(chunk),
           "--min-score", str(bit_score_cutoff),
           "--comp-based-stats", str(0),
           "--algo", str(0)]

    if sensitivity:
        cmd.append(f"--{sensitivity}")

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
        
    logger.info(f"[SEARCH_DONE] out={output_path}")
    return output_path, column_info_path



# ==========================
# search_filtering function
# ==========================
def filter_best_hits_by_score(search_output_path: str,
                     aln_score_cutoff: int):
    """
    Filters DIAMOND search output, keeping only the highest scoring hit per query above a minimum score threshold.

    Args:
        search_output_path (str): Path to DIAMOND search output file
        aln_score_cutoff (int): Minimum raw alignment score threshold to keep hits

    Returns:
        search_output_path (str): Path to filtered search output file
    """
    logger.info(f"[FILTER_START] Filtering search results keep_highest_score_per_query=True  aln_score_cutoff={aln_score_cutoff}")
    
    # Determine columns structure from the file content
    columns = ["qseqid", "sseqid", "pident", "qlen", "slen", "length", "mismatch", "gapopen", "qstart", "qend",
               "sstart", "send", "evalue", "bitscore", "score","qtitle"]
    
    best_hits = {}  # {qseqid: (score, full_row)}
    score_idx = columns.index('score')
    qseqid_idx = columns.index('qseqid')
    qtitle_idx = columns.index('qtitle')
    
    # Stream through file to find best hits
    with open(search_output_path, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            qseqid = fields[qseqid_idx]
            qtitle = fields[qtitle_idx]
            score = float(fields[score_idx])
            
            # Only keep hits with score >= threshold
            if score >= aln_score_cutoff:
                if qtitle not in best_hits or score > best_hits[qtitle][0]:
                    best_hits[qtitle] = (score, line)
    
    # ===========================
    # Write filtered results
    # ===========================
    temp_output = search_output_path + '.tmp'
    with open(temp_output, 'w') as f:
        for score, line in best_hits.values():
            f.write(line)
    
    # Replace original with filtered
    Path(temp_output).replace(search_output_path)
    
    logger.info(f"[FILTER_DONE] total_retained={len(best_hits)} out={search_output_path}")
    return search_output_path



# =============================
# parse_search_output function
# =============================
def split_search_output_by_mapping(blast_out_path: str,
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
    dict_seqid_group_name = {}
    for group_name, seq_ids in mapping_dict.items():
        for seq_id in seq_ids:
            if seq_id not in dict_seqid_group_name:
                dict_seqid_group_name[seq_id] = set()
            dict_seqid_group_name[seq_id].add(group_name)
    
    db_search_name = determine_dataset_name(blast_out_path, '.', 0, '').rsplit('_hits', 1)[0]
    
    logger.info(
        f"[SPLIT_START] source={blast_out_path} categories={len(mapping_dict)} "
        f"unique_keys_across_all_categories={len(dict_seqid_group_name):,} split_column={split_column}"
    )
    
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
    lines_with_insufficient_columns = 0
    unmatched_source_lines = 0
    
    try:
        with open(blast_out_path, 'r') as infile:
            for line in infile:
                total_lines += 1
                fields = line.strip().split('\t')

                if len(fields) <= split_column:
                    lines_with_insufficient_columns += 1
                    logger.warning(f"Line {total_lines} has insufficient columns for splitting: {line.strip()}")
                    continue

                seq_key = fields[split_column]

                if seq_key in dict_seqid_group_name.keys():
                    for group_name in dict_seqid_group_name[seq_key]:
                        file_handles[group_name].write(line)
                        line_counts[group_name] = line_counts.get(group_name, 0) + 1
                else:
                    unmatched_source_lines += 1

    except Exception as e:
        logger.error(f"Error while splitting BLAST output: {e}")
        raise RuntimeError(f"Failed to split BLAST output: {e}") from e
    
    # =======================
    # Close all file handles
    # =======================
    finally:
        for file_handle in file_handles.values():
            file_handle.close()
    
    logger.info(
        f"[SPLIT_DONE] source={blast_out_path} total_lines={total_lines:,} "
        f"insufficient_cols={lines_with_insufficient_columns:,} "
        f"source_lines_not_in_any_category={unmatched_source_lines:,}"
    )
    
    for name, count in line_counts.items():
        if count > 0:
            logger.info(f"[SPLIT_GROUP] group={name} hits={count:,} out={split_outputs[name]}")
        else:
            logger.info(f"[SPLIT_GROUP] group={name} hits=0")
    
    return split_outputs



# ====================================
# aastk extract_hit_seqs CLI FUNCTION
# ====================================
def extract_matched_reads(blast_tab: str,
                        query_path: str,
                        output_dir: str,
                        query_key_column: int=15,
                        force: bool = False):
    """
    Extracts read sequences that have hits in the DIAMOND BLAST tabular output.

    Args:
        blast_tab (str): Path to DIAMOND BLAST tabular output file
        query_path (str): Path to sequencing read file, can be gzipped
        output_dir (str): Directory to save output files
        query_key_column (int): Column index for query key (default: 15 = qtitle)
        force (bool): Whether to overwrite existing files

    Returns:
        hit_seqs_path (str): Path to the extracted hit sequences in gzipped FASTQ format
        id_file (str): Path to text file containing extracted IDs
        dataset_dict (dict): Dictionary with dataset name and list of query keys
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
    # Extract unique keys from dereplicated BLAST results and fill mapping dict for downstream use
    logger.info(f"[HIT_EXTRACT_START] source={blast_tab} query_key_column={query_key_column}")

    unique_keys = extract_unique_keys(blast_tab, query_key_column)

    dataset_name = determine_dataset_name(query_path, '.', 0, '')
    dataset_dict = {}
    dataset_dict[dataset_name] = list(unique_keys)

    # ======================================================
    # Write unique keys to file for seqkit grep
    # ======================================================
    id_file = ensure_path(output_dir, f"{protein_name}_matched_ids.txt", force=force)
    with open(id_file, 'w') as f:
        for key in unique_keys:
            f.write(f"{key}\n")

    # ======================================================
    # Extract and write matching sequences in gzipped FASTQ
    # ======================================================
    cmd = ["seqkit", "grep", "-n", "-f", id_file, query_path, "-o", out_fastq]

    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    
    except subprocess.CalledProcessError as e:
        logger.error(f"seqkit grep failed: {e.stderr}")
        raise RuntimeError(f"Failed to extract matching sequences: {e.stderr}") from e
    
    logger.info(
        f"[HIT_EXTRACT_DONE] keys={len(unique_keys)} source={blast_tab} ids={id_file} out={out_fastq}"
    )

    return out_fastq, id_file, dataset_dict




# ================================
# merge_extracted_hit_seqs function
# ================================
def merge_matched_reads(all_hit_seqs_paths: list,
                output_dir: str,
                force: bool = False):
    """
    Merges multiple FASTQ files containing hit sequences into a single FASTQ file.
    Deduplication is done by read name (not by sequence) to preserve dataset-specific
    IDs when different datasets contain identical read sequences.

    Args:
        all_hit_seqs_paths (list): List of paths to FASTQ files containing hit sequences
        output_dir (str): Directory to save merged FASTQ file
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
        
    # =======================================================
    # Merge and deduplicate all matched sequences (search 1)
    # =======================================================
    logger.info(f"Merging {len(all_hit_seqs_paths)} hit FASTQ files with name-based deduplication")

    cmd = ["seqkit", "rmdup", "-n", "-o", merged_hits_file, *all_hit_seqs_paths]

    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        logger.info(f"Successfully merged and name-deduplicated {len(all_hit_seqs_paths)} FASTQ files into {merged_hits_file}")
    
    except subprocess.CalledProcessError as e:
        logger.error(f"seqkit rmdup failed: {e.stderr}")
        raise RuntimeError(f"Failed to merge and deduplicate hit sequences: {e.stderr}") from e
    
    return merged_hits_file, infile_list




# ===============================
# bsr function
# ===============================
def compute_bsr(seed_search_tab: str,
        og_search_tab: str,
        column_info_path: str,
        output_dir: str,
        force: bool = False):
    """
    Computes BSR (Blast Score Ratio) using protein BLAST tab file and outgroup BLAST tab file.

    Args:
        seed_search_tab (str): Path to BLAST output file from Search 1 (seed db)
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
    protein_name = determine_dataset_name(seed_search_tab, '.', 0, '_hits')
    bsr_file = ensure_path(output_dir, f"{protein_name}_bsr.tsv", force=force)

    logger.info(
        f"[BSR_COMPUTE_START] target={protein_name} seed_hits={seed_search_tab} og_hits={og_search_tab}"
    )
    
    # ==========================================================
    # Load column info if provided and retrieve columns indexes
    # ==========================================================
    if column_info_path and Path(column_info_path).exists():
        try:
            with open(column_info_path, 'r') as f:
                column_info = json.load(f)
                
                # Get all required column indexes from .json file
                qseqid_col = column_info.get('qseqid', 0)
                qtitle_col = column_info.get('qtitle', 15)
                sseqid_col = column_info.get('sseqid', 1)
                pident_col = column_info.get('pident', 2)
                length_col = column_info.get('length', 5)
                score_column = column_info.get('score', 14)
                
                logger.info(f"[BSR_COLUMN_INFO] Column indexes: qseqid={qseqid_col}, qtitle={qtitle_col}, sseqid={sseqid_col}, "
                           f"pident={pident_col}, length={length_col}, score={score_column}")

        except Exception as e:
            logger.warning(f"Failed to load column info from {column_info_path}: {e}")
            # Fall back to default indexes
            qseqid_col = 0
            qtitle_col = 15
            sseqid_col = 1
            pident_col = 2
            length_col = 5
            score_column = 14
    else:
        # Use default column indexes if no .json file provided
        qseqid_col = 0
        qtitle_col = 15
        sseqid_col = 1
        pident_col = 2
        length_col = 5
        score_column = 14
        logger.info("[BSR_COLUMN_INFO] Using default column indexes (no column_info file provided)"
                   f" qseqid={qseqid_col}, qtitle={qtitle_col}, sseqid={sseqid_col}, "
                   f"pident={pident_col}, length={length_col}, score={score_column}")

    # =====================================
    # Load BLAST outputs and calculate BSR
    # =====================================
    # Read outgroup results into dict
    db_data = {}  # {qseqid: (qtitle, sseqid_db, pident_db, length_db, score_db)}
    total_seed_count = 0
    
    with open(seed_search_tab, 'r') as f:
        for line in f:
            total_seed_count += 1
            fields = line.strip().split('\t')
            qseqid = fields[qseqid_col]
            qtitle = fields[qtitle_col]
            sseqid_db = fields[sseqid_col]
            pident_db = float(fields[pident_col])
            length_db = float(fields[length_col])
            score_db = float(fields[score_column])
            db_data[qseqid] = (qtitle, sseqid_db, pident_db, length_db, score_db)
    
    # Stream through outgroup results and calculate BSR
    matched_count = 0
    total_og_count = 0
    
    with open(bsr_file, 'w') as out_f:
        # Write header
        out_f.write('qseqid\tqtitle\tsseqid_db\tpident_db\tlength_db\tscore_db\tsseqid_og\tpident_og\tlength_og\tscore_og\tBSR\n')
        
        with open(og_search_tab, 'r') as og_f:
            for line in og_f:
                fields = line.strip().split('\t')
                qseqid = fields[qseqid_col]
                total_og_count += 1
                
                # Only process if we have database data for this query
                if qseqid in db_data:
                    sseqid_og = fields[sseqid_col]
                    pident_og = float(fields[pident_col])
                    length_og = float(fields[length_col])
                    score_og = float(fields[score_column])
                    
                    qtitle, sseqid_db, pident_db, length_db, score_db = db_data[qseqid]
                    
                    # Calculate BSR
                    bsr_value = score_db / score_og if score_og > 0 else 0.0
                    
                    # Write result
                    out_f.write(f"{qseqid}\t{qtitle}\t{sseqid_db}\t{pident_db}\t{length_db}\t{score_db}\t{sseqid_og}\t{pident_og}\t{length_og}\t{score_og}\t{bsr_value}\n")
                    matched_count += 1
    
    logger.info(
        f"[BSR_COMPUTE_DONE] target={protein_name} matched={matched_count:,} total_seed={total_seed_count:,} total_og={total_og_count:,} out={bsr_file}"
    )
    
    return bsr_file




# ===============================
# aastk rasr_plot CLI function
# ===============================
def rasr_plot(bsr_file: str,
                output_dir: str,
                bsr_cutoff: float = 0.9,
                dbmin: int = None,
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

    logger.info(f"[BSR_PLOT_START] Creating BSR scatter plot for {protein_name}")

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
        logger.info(f"[BSR_PLOT_DONE] BSR scatter plot saved to {out_graph}")

        return out_graph

    except Exception as e:
        logger.error(f"Error creating BSR plot: {e}")
        raise RuntimeError(f"Failed to create BSR plot: {e}") from e




# ===============================
# aastk rasr_select CLI function
# ===============================
def rasr_select(dbmin: float,
                bsr_cutoff: float,
                matched_fastq: str,
                bsr_file: str,
                output_dir: str,
                force: bool = False):
    """
    Subsets RASR results based on BSR and database score thresholds.

    Args:
        dbmin (float): Minimum database score threshold for selection
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

    # Handle None values from cli.py
    if bsr_cutoff is None:
        bsr_cutoff = 0.9
    if dbmin is None:
        dbmin = 110

    logger.info(f"[SELECT_START] Selecting RASR hits with database score >= {dbmin} and BSR >= {bsr_cutoff}")

    # ===============================
    # Output file path setup
    # ===============================
    protein_name = determine_dataset_name(bsr_file, '.', 0, '_bsr')
    out_fastq = ensure_path(output_dir, f"{protein_name}_selected.fastq", force=force)
    id_file = ensure_path(output_dir, f"{protein_name}_selected_ids.txt", force=force)
    out_tsv = ensure_path(output_dir, f"{protein_name}_selected.tsv", force=force)
    
    # ===================================
    # Extract reads IDs based on cutoffs
    # ===================================
    # Subset the bsr table by the cutoffs, save it and write selected IDs to file
    bsr_df = pd.read_csv(bsr_file, sep='\t', header=0)

    selected_tsv = bsr_df[(bsr_df['score_db'] >= dbmin) & (bsr_df['BSR'] >= bsr_cutoff)]
    selected_tsv.to_csv(out_tsv, sep='\t', index=False)

    selected_ids = selected_tsv['qseqid'].dropna().astype(str).unique()

    with open(id_file, 'w') as f:
        for seq_id in selected_ids:
            f.write(f"{seq_id}\n")

    logger.info(f"[SELECT_DONE] total_ids={len(selected_ids)} ids_out={id_file}")

    # ===============================
    # Extract selected sequences
    # ===============================
    try:
        cmd = ["seqkit", "grep", "-f", id_file, matched_fastq, "-o", out_fastq]

        subprocess.run(cmd, check=True, capture_output=True, text=True)
        logger.info(f"[SELECT_DONE] fastq_out={out_fastq}")

    
    except Exception as e:
        logger.error(f"seqkit grep failed: {e.stderr}")
        raise RuntimeError(f"Failed to extract selected sequences: {e.stderr}") from e

    return out_fastq, id_file, out_tsv




# ===============================
# aastk rasr WORKFLOW
# ===============================
def rasr(query: str,
            seed: str,
            query_dir: str,
            seed_dir: str,
            outgrp_db: str,
            output_dir: str,
            sensitivity: str,
            block: int,
            chunk: int,
            bit_score_cutoff: int,
            aln_score_cutoff: int,
            bsr_cutoff: float,
            dbmin: int,
            threads: int = 1,
            force: bool = False,
            keep: bool = False):
    """
    Run the RASR workflow.

    Pipeline overview:
        1. Validate inputs and build execution plan.
        2. Optionally merge seed FASTA files (batch seed mode).
        3. Build DIAMOND seed database.
        4. For each query dataset:
            - search seed DB
            - filter best hit per query (alignment score cutoff)
            - optionally split hits per gene (batch seed mode)
            - extract matched reads per dataset x gene
        5. Optionally merge matched reads across datasets/genes. (batch modes)
        6. Search outgroup database with merged/single matched reads.
        7. Optionally split outgroup hits per dataset (batch query mode).
        8. For each dataset x gene pair:
            - compute BSR
            - generate BSR plot
            - select reads passing dbmin + bsr_cutoff

    Args:
        query (str): Path to a single query FASTQ file (mutually exclusive with query_dir).
        seed (str): Path to a single seed FASTA file (mutually exclusive with seed_dir).
        query_dir (str): Path to a directory of query FASTQ files.
        seed_dir (str): Path to a directory of seed FASTA files.
        outgrp_db (str): Path to outgroup DIAMOND database.
        output_dir (str): Directory for output files.
        sensitivity (str): DIAMOND sensitivity preset (for example: fast, sensitive).
        block (int): DIAMOND block size parameter.
        chunk (int): DIAMOND chunk size parameter.
        bit_score_cutoff (int): DIAMOND minimum score cutoff (default: 10).
        aln_score_cutoff (int): Post-search alignment score cutoff (default: 50).
        threads (int): Number of CPU threads.
        force (bool): Overwrite existing files.
        keep (bool): Keep intermediate files.
        bsr_cutoff (float): Minimum BSR threshold for final selection (default: 0.9).
        dbmin (int): Minimum database score threshold for final selection (default: 110).

    Returns:
        dict: Mapping {(dataset_name, gene_name): result_info} where result_info includes
            key output file paths (search output, matched reads, BSR table, plot, selected reads).
    """

    logger.info("[RUN_START] workflow=rasr")
    logger.info("[INITILISATION_START] validating inputs, building execution plan and diamond seed database")

    # =====================
    # Build execution plan
    # =====================
    logger.info("[PLAN_START] builder=build_execution_plan")
    plan = build_execution_plan(query, query_dir, seed, seed_dir)
    
    # =======================================
    # List inputs and plan mode of execution
    # =======================================
    logger.info("[LISTING_INPUTS_START]")

    query_files = list_inputs(plan["query_info"])
    db_files = list_inputs(plan["seed_info"])
    
    if not db_files:
        raise ValueError("No seed database files found")
    if not query_files:
        raise ValueError("No query files found")
    
    # =====================================
    # Output directory and variables setup
    # =====================================
    logger.info(
        f"[RUN_CONFIG] query_datasets={len(query_files)} seed_dbs={len(db_files)} "
        f"query_mode={plan['query_info']['mode']} seed_mode={plan['seed_info']['mode']} out={output_dir}"
    )

    intermediate_results = {}
    results_by_ds_gene = {}
    db_seqid_dict = None

    try:
        # ==================================
        # Gene of interest database merging
        # ==================================
        if plan['steps']['merge_seed_databases']:
            logger.info(f"[MERGE_SEED_DBS_START] seed_files={len(db_files)}")

            merged_seed_db_path, db_seqid_dict = merge_seed_databases(db_files, output_dir, force=force)
            db_input_for_build = merged_seed_db_path

            intermediate_results['merged_seed_db_fasta'] = merged_seed_db_path

        else:
            logger.info("[MERGE_SEED_DBS_SKIP] reason=single seed input")

            db_input_for_build = str(db_files[0])

        # ===================================
        # Gene of interest database building
        # ===================================
        logger.info(f"[BUILD_DMND_DB_START] source={db_input_for_build}")

        db_path = pasr_build(db_input_for_build, threads, output_dir, force=force)

        logger.info(f"[BUILD_DMND_DB_DONE] out={db_path}")

        intermediate_results['merged_seed_db_diamond'] = db_path

        # ===================================
        # Track outputs per dataset and gene
        # ===================================
        all_hit_seqs_paths = []  # Accumulate all _matched.fastq.gz file paths across datasets and genes
        dict_dataset_seqids = {}  # {dataset_name: [qseqids]}

        # ==================================================
        # Step 1: Search each dataset against gene database
        # ==================================================     
        logger.info(f"[STEP_1_START] search_seed_db datasets={len(query_files)} genes={len(db_files)}")

        dataset_output_dirs = set()
        for query_file in query_files:
            dataset_name = determine_dataset_name(query_file, '.', 0, '')
            dataset_output_dir = ensure_path(output_dir, dataset_name, force=force)
            dataset_output_dirs.add(dataset_output_dir)
            
            logger.info(f"[DATASET_START] workflow_step=1 dataset={dataset_name}")

            # ======================
            # Search gene database
            # ======================
            logger.info(f"[SEARCH_1_START] dataset={dataset_name} query={query_file} db={db_path}")
            search_output, column_info_path = search(db_path, str(query_file), threads, dataset_output_dir, sensitivity, bit_score_cutoff=bit_score_cutoff, max_target_seqs=1, block=block, chunk=chunk, force=force)
            
            # =========================
            # Filter search results
            # =========================
            search_output = filter_best_hits_by_score(search_output, aln_score_cutoff=aln_score_cutoff)

            if 'search_output' not in intermediate_results:
                intermediate_results['search_output'] = []
            intermediate_results['search_output'].append(search_output)

            # ==============================
            # Split search results per gene
            # ==============================
            if plan['steps']['parse_search_1_per_db']:
                if db_seqid_dict is None:
                    raise RuntimeError("mapping dict database:[seqid] is required for splitting search outputs when seed_is_multi is True but was not initialized")

                db_split_outputs = split_search_output_by_mapping(search_output, mapping_dict=db_seqid_dict, output_dir=dataset_output_dir, split_column=1, force=force)

                # Move split outputs to gene-specific subdirectories
                for gene_name, split_output in list(db_split_outputs.items()):
                    gene_output_dir = ensure_path(dataset_output_dir, gene_name, force=force)
                    dest_path = ensure_path(gene_output_dir, Path(split_output).name, force=force)
                    Path(split_output).replace(dest_path)
                    db_split_outputs[gene_name] = dest_path

            else:
                single_gene_name = determine_dataset_name(str(db_files[0]), '.', 0, '')
                db_split_outputs = {single_gene_name: search_output}


            if 'db_split_outputs' not in intermediate_results:
                intermediate_results['db_split_outputs'] = []
            intermediate_results['db_split_outputs'].append(db_split_outputs)            

            # ===============================
            # Extract hit sequences per gene
            # ===============================
            for gene_name, split_output in db_split_outputs.items():
                gene_output_dir = ensure_path(dataset_output_dir, gene_name, force=force) 

                hit_seqs_path, id_file, dataset_dict = extract_matched_reads(split_output, str(query_file), 
                                                                        output_dir=gene_output_dir,
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
                if dataset_name not in dict_dataset_seqids:
                    dict_dataset_seqids[dataset_name] = []
                dict_dataset_seqids[dataset_name].extend(dataset_dict[dataset_name])

            logger.info(f"[DATASET_DONE] step=1 dataset={dataset_name} genes_processed={len(db_split_outputs)}")

        if len(dataset_output_dirs) != len(query_files):
            raise RuntimeError(
                "Dataset output directory count does not match input datasets: "
                f"{len(dataset_output_dirs)} dirs for {len(query_files)} datasets"
            )

        logger.info(f"[STEP_1_DONE] datasets_processed={len(dataset_output_dirs)}")

        intermediate_results['column_info_path'] = column_info_path

        # ============================================================================
        # Step 2: Merge hits from all datasets and genes and search outgroup database
        # ============================================================================
        logger.info("[STEP_2_START] search_outgroup")

        do_merge_hits = plan['steps']['merge_extracted_hit_seqs'] and len(all_hit_seqs_paths) > 1

        if do_merge_hits:
            logger.info(f"[HIT_MERGE_START] files={len(all_hit_seqs_paths)}")
            merged_hits_file, infile_list = merge_matched_reads(all_hit_seqs_paths, output_dir, force=force)
            intermediate_results['merged_hits_file'] = merged_hits_file
            intermediate_results['infile_list'] = infile_list
        else:
            logger.info("[HIT_MERGE_SKIP] reason=single_input_fastq")
            merged_hits_file = all_hit_seqs_paths[0]

        # =================== Search outgroup database with merged hits ===================

        logger.info(f"[SEARCH_2_START] query={merged_hits_file}")
        og_search_output, og_column_info_path = search(outgrp_db, merged_hits_file, threads, output_dir, sensitivity, bit_score_cutoff=bit_score_cutoff, max_target_seqs=25, block=block, chunk=chunk, force=force)

        intermediate_results['og_search_output'] = og_search_output
        intermediate_results['og_column_info_path'] = og_column_info_path

        # =================== Filter outgroup search results ===================

        og_search_output_filtered = filter_best_hits_by_score(og_search_output, aln_score_cutoff=aln_score_cutoff)

        intermediate_results['og_search_output_filtered'] = og_search_output_filtered

        # =================== Split outgroup search results per dataset ===================

        og_split_outputs = {}
        if plan['steps']['parse_search_2_per_query']:
            logger.info(f"[SPLIT_2_START] datasets={len(dict_dataset_seqids)} split_column=15")

            og_split_outputs = split_search_output_by_mapping(og_search_output_filtered, mapping_dict=dict_dataset_seqids, output_dir=output_dir, split_column=15, force=force)

            # Move split outputs to dataset-specific subdirectories
            for dataset_name, split_output in list(og_split_outputs.items()):
                dataset_output_dir = ensure_path(output_dir, dataset_name, force=force)
                dest_path = ensure_path(dataset_output_dir, Path(split_output).name, force=force)
                Path(split_output).replace(dest_path)
                og_split_outputs[dataset_name] = dest_path
                
        else:
            logger.info("[SPLIT_2_SKIP] reason=single_dataset")

            only_dataset = next(iter(dict_dataset_seqids.keys()))
            og_split_outputs[only_dataset] = og_search_output_filtered

        intermediate_results['og_split_outputs'] = og_split_outputs
        logger.info(f"[STEP_2_DONE] outgroup_splits={len(og_split_outputs)}")

        # ================================================================
        # Step 3: Calculate BSR and generate plots for each dataset x gene
        # ================================================================
        logger.info(f"[STEP_3_START] pairs={len(results_by_ds_gene)}")

        for (dataset_name, gene_name), result_info in results_by_ds_gene.items():
            
            gene_output_dir = result_info['output_dir']
            db_split_output = result_info['db_search_output']
            og_split_output = og_split_outputs.get(dataset_name)
            
            if db_split_output and og_split_output:
                logger.info(f"[PAIR_START] dataset={dataset_name} gene={gene_name}")

                # =================== Calculate BSR values ===================
                bsr_output = compute_bsr(db_split_output, og_split_output, column_info_path, output_dir=gene_output_dir, force=force)

                result_info['bsr_output'] = bsr_output

                # =================== Plot BSR values ===================
                plot_output = rasr_plot(bsr_output, output_dir=gene_output_dir, bsr_cutoff=bsr_cutoff, dbmin=dbmin, force=force)

                result_info['bsr_plot'] = plot_output
          
                # =================== Select hits above BSR and dbmin thresholds ===================
                hit_seqs_path = result_info['matched_fastq_path']
                selected_fastq, selected_ids_file, selected_tsv_file = rasr_select(dbmin, bsr_cutoff, str(hit_seqs_path), bsr_output, output_dir=gene_output_dir, force=force)

                result_info['selected_fastq'] = selected_fastq
                result_info['selected_tsv_file'] = selected_tsv_file
                intermediate_results['selected_ids_file'] = selected_ids_file
                logger.info(f"[PAIR_DONE] dataset={dataset_name} gene={gene_name} selected_ids={selected_ids_file}")
            else:
                logger.warning(f"[PAIR_SKIP] dataset={dataset_name} gene={gene_name} reason=missing_split_outputs")

        logger.info("[STEP_3_DONE] pair_processing_complete=true")
        logger.info("[RUN_DONE] workflow=rasr_multiple status=success")
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
                    
