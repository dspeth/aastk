#!/usr/bin/env python3

from .util import extract_unique_keys, determine_file_type, write_fa_matches, write_fq_matches

import subprocess
import os
import matplotlib.pyplot as plt
import pandas as pd
from io import StringIO

def build_protein_db(protein_name: str, seed_fasta: str, threads: int, db_dir: str):
    """
        Builds a DIAMOND protein database from a seed FASTA file.

        Args:
            protein_name (str): Name of the protein for the database.
            seed_fasta (str): Path to the FASTA file containing seed sequences.
            threads (int): Number of threads to use.
            db_dir (str): Directory where the database should be stored. (Default: current working directory)

        Returns:
            db_path: Path to DIAMOND protein database.

        Raises:
            RuntimeError: If the DIAMOND database creation fails.
    """
    # check for db_dir
    if db_dir is None:
        db_dir = '.'

    # configure target directory and ensure target directory exists
    os.makedirs(db_dir, exist_ok=True)
    db_path = os.path.join(db_dir, f"{protein_name}_seed_db")

    #try running the subprocess for the Diamond makedb command
    try:
        subprocess.run(
            ["diamond", "makedb", "--in", seed_fasta, "-d", db_path, "-p", str(threads)],
            check=True
        )
    except subprocess.CalledProcessError as e:
        print(f"Error in building the DIAMOND database: {e}")
        raise

    return db_path

def search_protein_db(db_path: str, query_path: str, protein_name: str, threads: int, output_dir: str, sensitivity: str, block: int, chunk: int):
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

    Returns:
        output_path: Path to tabular BLAST output file.
    """
    # check for output_dir
    if output_dir is None:
        output_dir = '.'

    # make sure target_dir exists and define results location and file name
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, f"{protein_name}_hits.txt")

    # define the output columns of interest
    columns = ["qseqid", "sseqid", "pident", "qlen", "slen", "length", "mismatch", "gapopen", "qstart", "qend",
               "sstart", "send", "evalue", "bitscore", "score"]

    # check for sensitivity, if None set to default --fast
    if sensitivity is None:
        sensitivity = '--fast'
    else:
        sensitivity = '--' + str(sensitivity)

    # check for block size
    if block is None:
        block = 6

    # check for chunk number
    if chunk is None:
        chunk = 2

    try:
        # run diamond blastp
        subprocess.run(
            ["diamond", "blastp", "-d", db_path, "-q", query_path, "-p", str(threads), "-o",
             output_path, "-k", str(1), "--matrix", "blosum45", "--masking", str(0), "--outfmt", str(6), *columns, "-b", str(block), "-c", str(chunk),
             "--min-score", str(50), "--comp-based-stats", str(0), sensitivity],
            check=True
        )
    except subprocess.CalledProcessError as e:
        print(f"Error in DIAMOND blast p search: {e}")
        raise

    return output_path

def extract_matching_sequences(blast_tab: str, query_path: str, output_dir: str, key_column: int = 0):
    """
    Extracts reads that have BLAST/DIAMOND hits against a custom database.

    Args:
        blast_tab: Tabular BLAST/DIAMOND output file.
        query_path: Fasta or fastq file containing sequencing reads used as BLAST/DIAMOND queries.
        output_dir (str): Directory where extracted sequences should be stored.
        key_column: Column index in the BLAST tab file to pull unique IDs from (default is 0).
    Returns:
        out_fasta: Path to output FASTA file.
    """
    # check for output_dir
    if output_dir is None:
        output_dir = '.'

    # make sure target directory exists and define path to output FASTA
    os.makedirs(output_dir, exist_ok=True)
    out_fasta = os.path.join(output_dir, "matched_sequences.fasta")

    # Extract unique keys (query IDs) from the specified column of the BLAST tab file
    matching_ids = extract_unique_keys(blast_tab, key_column)

    # Determine file type (fasta or fastq)
    file_type = determine_file_type(query_path)

    # Open the output file for writing
    with open(out_fasta, "w") as out:
        # Use the appropriate generator based on the file type
        if file_type == "fasta":
            print("File type: FASTA")
            for header, sequence in write_fa_matches(query_path, matching_ids):
                out.write(f"{header}\n{sequence}\n")
        elif file_type == "fastq":
            print("File type: FASTQ")
            for header, sequence in write_fq_matches(query_path, matching_ids):
                out.write(f"{header}\n{sequence}\n")

    return out_fasta

def calculate_max_scores(extracted: str, matrix: str, output_dir: str):
    """
    Calculates max scores for sequences using a BLOSUM matrix.

    Args:
        extracted (str): Path to the extracted FASTA file.
        matrix (str): BLOSUM matrix name ('BLOSUM45' or 'BLOSUM62').
        output_dir (str): Directory where the output file should be stored.

    Returns:
        max_scores: Dictionary of protein headers and their max scores.
    """
    # check for output_dir
    if output_dir is None:
        output_dir = '.'

    # check if target_dir exists and define path to output file
    os.makedirs(output_dir, exist_ok=True)
    out_file = os.path.join(output_dir, "max_scores.txt")

    # check if valid matrix was chosen
    if matrix not in ["BLOSUM45", "BLOSUM62"]:
        raise ValueError("Matrix must either be 'BLOSUM45' or 'BLOSUM62'!")

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

    # read the fasta file containing matched sequences as obtained via pasr extract
    sequences = {}
    current_header = None

    with open(extracted, 'r') as file:
        for line in file:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                # get sequence ID and store corresponding sequences in the sequences dictionary
                current_header = line.replace('>','')
                sequences[current_header] = ""
            else:
                if current_header:
                    sequences[current_header] += line

        # calculate the max scores
        max_scores = {}
        for header, sequence in sequences.items():
            score = 0
            for amino_acid in sequence:
                if amino_acid in blosum_diagonals[matrix]:
                    score += blosum_diagonals[matrix][amino_acid]
            max_scores[header] = score

        # write the results of the step to specified output file
        with open(out_file, 'w') as out:
            out.write("Protein_id\tmax_score\n")
            for header, score in max_scores.items():
                out.write(f"{header}\t{score}\n")

    return out_file

def blast_score_ratio(blast_tab: str, max_scores_path: str, output_dir: str, key_column: int = 0):
    """
    Computes BSR (Blast Score Ratio) using a BLAST tab file and max scores from a TSV.

    Args:
        blast_tab (str): Path to DIAMOND/BLAST output file (must include 'score' column).
        max_scores_path (str): Path to TSV file with max scores (headers: Protein_id, max_score).
        output_dir (str): Directory to save the BSR results.
        key_column (int): Column index in blast_tab to use for matching. (Default: 0 for 'qseqid')

    Returns:
        bsr_output (str): Path to the output file with BSR values.
    """
    # check for output_dir
    if output_dir is None:
        output_dir = '.'

    # ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    bsr_file = os.path.join(output_dir, "blast_score_ratios.txt")

    # parse the max_scores.tsv file
    max_scores = {}
    with open(max_scores_path) as tsv:
        header = tsv.readline()
        for line in tsv:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                max_scores[parts[0]] = float(parts[1])

    with open(blast_tab) as infile, open(bsr_file, 'w') as out:
        out.write("qseqid\tsseqid\tscore\tmax_score\tBSR\n")

        for line in infile:
            parts = line.strip().split('\t')
            key = parts[key_column]
            try:
                raw_score = float(parts[14])
                max_score = max_scores[key]
                bsr = raw_score / max_score
                out.write(f"{parts[0]}\t{parts[1]}\t{raw_score:.1f}\t{max_score:.1f}\t{bsr:.4f}\n")
            except KeyError:
                print(f"Warning: No max score found for {key}")
            except ValueError:
                print(f"Warning: Couldn't convert score to float for line: {line.strip()}")

    return bsr_file

def plot_bsr(bsr_file: str, output_dir: str):
    # check for output_dir
    if output_dir is None:
        output_dir = '.'

    # check if target_dir exists and define path to output file
    os.makedirs(output_dir, exist_ok=True)
    out_graph = os.path.join(output_dir, 'bsr_plot.png')

    bsr_df = pd.read_csv(bsr_file, sep='\t', header=0)

    scatter = plt.scatter(
        bsr_df['max_score'],
        bsr_df['score'],
        c=bsr_df['BSR'],
        cmap='viridis',
        edgecolor='k',
        alpha=0.75
    )

    plt.xlabel('Calculated maximum score')
    plt.ylabel('Alignment score to seed set')
    plt.xlim(0, 1.5 * bsr_df['max_score'].max())
    plt.axvline(bsr_df['max_score'].max(), color='red', linestyle='--', linewidth=1.2, label='Max score')
    plt.title('Blast Score Ratio (BSR) Scatter Plot')
    cbar = plt.colorbar(scatter)
    cbar.set_label('% seq identity')

    plt.tight_layout()
    plt.savefig(out_graph)
    plt.close()

def metadata(matched_seqs: str, selfmin: int, selfmax: int, metadata_file: str, dbmin: int = None,
             bsr: float = None):
    result = subprocess.run(
        ["seqkit", "stats", matched_seqs], capture_output=True, text=True
    )

    df = pd.read_csv(StringIO(result.stdout), sep="\t")

    df["selfmin"] = selfmin
    df["selfmax"] = selfmax
    df["dbmin"] = dbmin
    df["bsr_cutoff"] = bsr

    df.to_csv(metadata_file, sep="\t", index=False)

def pasr(db_dir: str, protein_name: str, seed_fasta: str, query_fasta: str, matrix_name: str,
         output_dir: str, sensitivity: str, block: int, chunk: int, key_column: int = 0, threads: int = 1,):
    """
    PASR workflow with configurable output directory.

    Args:
        db_dir (str): Path to DIAMOND database directory
        protein_name (str): Protein name.
        seed_fasta (str): Path to seed FASTA.
        query_fasta (str): Path to query FASTA.
        matrix_name (str): BLOSUM matrix ('BLOSUM45' or 'BLOSUM62').
        output_dir (str): Output directory for (default: current directory).
        sensitivity (str): Choose sensitivity of diamond blastp search (Default: --fast)
        block (int): Choose diamond blastp sequence block size in billions of letters. (Default: 6)
        chunk (int): Choose number of chunks for diamond blastp index processing. (Default: 2)
        key_column: Column index in the BLAST tab file to pull unique IDs from (default is 0).
        threads (int): Number of threads (default: 1).
    """
    # check for output_dir
    if output_dir is None:
        output_dir = '.'

    #get full path to target dir
    output_dir = os.path.abspath(output_dir)
    print(f"Running PASR workflow. Output directory: {output_dir}")

    db_path = build_protein_db(protein_name, seed_fasta, threads, db_dir)
    search_output = search_protein_db(db_path, query_fasta, protein_name, threads, output_dir, sensitivity, block, chunk)
    matched_fasta = extract_matching_sequences(search_output, query_fasta, output_dir, key_column)
    max_scores = calculate_max_scores(matched_fasta, matrix_name, output_dir)
    bsr_file = blast_score_ratio(search_output, max_scores, output_dir, key_column)
    return plot_bsr(bsr_file, output_dir)

