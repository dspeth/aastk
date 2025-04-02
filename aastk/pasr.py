#!/usr/bin/env python3

from .util import extract_unique_keys, determine_file_type

import subprocess
import os

def build_protein_db(protein_name: str, seed_fasta: str, threads: int):
    """
        Builds a DIAMOND protein database from a seed FASTA file.

        Args:
            protein_name (str): Name of the protein for the database.
            seed_fasta (str): Path to the FASTA file containing seed sequences.
            threads (int): Number of threads to use.

        Raises:
            RuntimeError: If the DIAMOND database creation fails.
    """
    # name the DB in accordance with pre-existing naming convention
    db_dir = f"{protein_name}_protein_db"
    db_name = f"{db_dir}/{protein_name}_seed_db"

    # ensure the directory exists
    os.makedirs(db_dir, exist_ok=True)

    #try running the subprocess for the Diamond makedb command
    try:
        subprocess.run(
            ["diamond", "makedb", "--in", seed_fasta, "-d", db_name, "-p", str(threads)],
            check=True
        )
    except subprocess.CalledProcessError as e:
        print(f"Error in building the DIAMOND database: {e}")
        raise

    return db_name

def search_protein_db(db_path: str, query_path: str, threads: int, protein_name: str):
    """
    Searches a DIAMOND reference database for homologous sequences.

    Args:
        db_path (str): Path to the DIAMOND reference database.
        query_path (str): Path to the query FASTA file containing sequences to be searched.
        threads (int): Number of CPU threads to use for the search.
        protein_name (str): Name of the target protein (used for output file naming or filtering).
    Returns:
        None: The function performs the search and outputs results to a file.
    """
    # define output path
    output_dir = f"{protein_name}_protein_db"
    output_path = f"{output_dir}/{protein_name}_hits"

    # ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    try:
        # run diamond blastp
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

def extract_matching_sequences(blast_tab: str, seq_file: str, out_fasta: str, key_column: int = 0):
    """
    Extracts reads that have BLAST/DIAMOND hits against a custom database.

    Args:
    - blast_tab: Tabular BLAST/DIAMOND output file.
    - read_file: Fasta or fastq file containing sequencing reads used as BLAST/DIAMOND queries.
    - out_fasta: Output file to store matched sequences.
    - key_column: Column index in the BLAST tab file to pull unique IDs from (default is 0).
    """
    # Extract unique keys (query IDs) from the specified column of the BLAST tab file
    matching_ids = extract_unique_keys(blast_tab, key_column)

    # Determine file type (fasta or fastq)
    file_type = determine_file_type(seq_file)


    # Open the output file for writing
    with open(out_fasta, "w") as out:
        # Use the appropriate generator based on the file type
        if file_type == "fasta":
            print("File type: FASTA")
            for header, sequence in write_fa_matches(seq_file, matching_ids):
                out.write(f">{header}\n{sequence}\n")
        elif file_type == "fastq":
            print("File type: FASTQ")
            for header, sequence in write_fq_matches(seq_file, matching_ids):
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

def calculate_max_scores(fasta_file: str, matrix_name: str, out_file: str):
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