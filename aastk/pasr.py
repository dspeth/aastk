#!/usr/bin/env python3

from .util import extract_unique_keys, determine_file_type

import subprocess
import os

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

def search_protein_db(db_path: str, query_path: str, protein_name: str, threads: int, output_dir: str):
    """
    Searches a DIAMOND reference database for homologous sequences.

    Args:
        db_path (str): Path to the DIAMOND reference database.
        query_path (str): Path to the query FASTA file containing sequences to be searched.
        threads (int): Number of CPU threads to use for the search.
        protein_name (str): Name of the target protein (used for output file naming or filtering).
        output_dir (str): Directory where results should be stored. (Default: current working directory)

    Returns:
        output_path: Path to output file.
    """
    # check for output_dir
    if output_dir is None:
        output_dir = '.'

    # make sure target_dir exists and define results location and file name
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, f"{protein_name}_hits.txt")

    try:
        # run diamond blastp
        # improvements: column names as single variable [expand list]
        # -b and -c should be configurable by user, select reasonable defaults (leave them as they are for now)
        # make sensitivity configurable
        subprocess.run(
            ["diamond", "blastp", "-d", db_path, "-q", query_path, "-p", str(threads), "-o",
             output_path, "-k", str(1), "--matrix", "blosum45", "--masking", str(0), "--outfmt", str(6),
             "qseqid", "sseqid", "pident", "qlen", "slen", "length", "mismatch", "gapopen", "qstart",
             "qend", "sstart", "send", "evalue", "bitscore", "score", "-b", str(6), "-c", str(2),
             "--min-score", str(50), "--comp-based-stats", str(0), "--sensitive"],
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
                out.write(f">{header}\n{sequence}\n")
        elif file_type == "fastq":
            print("File type: FASTQ")
            for header, sequence in write_fq_matches(query_path, matching_ids):
                out.write(f">{header}\n{sequence}\n")

    return out_fasta

def write_fa_matches(seq_file, ids):
    """
    Generator function to process FASTA file and yield matching sequences in fasta format.

    Args:
    - seq_file: Path to the fasta file containing the sequences to search.
    - ids: set of ids to retrieve matches for.

    Yields:
    - Header and sequence of matching sequences.
    """
    matching = False
    sequence = ""

    ## parser can also be in util.py
    with open(seq_file, 'r') as sf:
        for line in sf:
            line = line.strip()
            if line.startswith(">"):
                if matching:
                    yield (header, sequence)
                sequence = ""
                seq_id = line.split()[0][1:]  # Get the query ID without '>'
                if seq_id in ids:
                    matching = True
                    header = line
                else:
                    matching = False
            elif matching:
                sequence += line

        if matching:
            yield (header, sequence)


def write_fq_matches(seq_file, ids):
    """
    Generator function to process FASTQ file and yield matching sequences in fasta format.

    Args:
    - seq_file: Path to the fastq file containing the sequences to search.
    - ids: set of ids to retrieve matches for.

    Yields:
    - Header and sequence of matching sequences (converted to fasta format).
    """
    matching = False
    line_count = 0
    sequence = ""

    with open(seq_file, 'r') as sf:
        for line in sf:
            line = line.strip()
            line_count += 1

            if line_count == 1:
                seq_id = line.split()[0][1:]  # Get the query ID without '@'
                matching = seq_id in ids
                if matching:
                    header = seq_id  # Store the fastq ID to convert to fasta format
            elif line_count == 2 and matching:
                sequence = line  # Store the sequence for matching read
            elif line_count == 4:
                line_count = 0  # Reset after each fastq record
                if matching:
                    yield (header, sequence)

def calculate_max_scores(fasta_file: str, matrix_name: str, output_dir: str):
    """
    Calculates max scores for sequences using a BLOSUM matrix.

    Args:
        fasta_file (str): Path to the extracted FASTA file.
        matrix_name (str): BLOSUM matrix name ('BLOSUM45' or 'BLOSUM62').
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
    if matrix_name not in ["BLOSUM45", "BLOSUM62"]:
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

    with open(fasta_file, 'r') as file:
        for line in file:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                # get sequence ID and store corresponding sequences in the sequences dictionary
                # get rid of ">"
                current_header = line.split()[0][1:]
                sequences[current_header] = ""
            else:
                if current_header:
                    sequences[current_header] += line

        # calculate the max scores
        max_scores = {}
        for header, sequence in sequences.items():
            score = 0
            for amino_acid in sequence:
                if amino_acid in blosum_diagonals[matrix_name]:
                    score += blosum_diagonals[matrix_name][amino_acid]
            max_scores[header] = score

        # write the results of the step to specified output file
        with open(out_file, 'w') as out:
            out.write("Protein_id\tmax_score\n")
            for header, score in max_scores.items():
                out.write(f"{header}\t{score}\n")

    return max_scores


def pasr(db_dir: str, protein_name: str, seed_fasta: str, query_fasta: str, matrix_name: str,
         output_dir: str, key_column: int = 0, threads: int = 1,):
    """
    PASR workflow with configurable output directory.

    Args:
        db_dir (str): Path to DIAMOND database directory
        protein_name (str): Protein name.
        seed_fasta (str): Path to seed FASTA.
        query_fasta (str): Path to query FASTA.
        matrix_name (str): BLOSUM matrix ('BLOSUM45' or 'BLOSUM62').
        output_dir (str): Output directory for (default: current directory).
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
    search_output = search_protein_db(db_path, query_fasta, protein_name, threads, output_dir)
    matched_fasta = extract_matching_sequences(search_output, query_fasta, output_dir, key_column)
    return calculate_max_scores(matched_fasta, matrix_name, output_dir)