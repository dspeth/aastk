#!/usr/bin/env python3

def extract_unique_keys(file_path, column_index=0):
    """
    Extracts and deduplicates keys from a specified column in a tab-delimited file.
    Args:
    - file_path: Path to the input file.
    - column_index: The index of the column to extract unique keys from (0-based).
    Returns:
    - A list of unique keys from the specified column.
    """
    unique_keys = set()
    with open(file_path, 'r') as file:
        for line in file:
            key = line.split('\t')[column_index]
            unique_keys.add(key.strip())
    return unique_keys

def determine_file_type(file_path):
    """
    Determines the file type (fasta or fastq) based on the first character in the file.
    Args:
    - file_path: Path to the file to check.
    Returns:
    - A string, either "fasta" or "fastq", indicating the file type.
    Raises:
    - ValueError: If the file type cannot be determined.
    """
    with open(file_path, 'r') as file:
        first_char = file.read(1)
        if first_char == '>':
            return "fasta"
        elif first_char == '@':
            return "fastq"
        else:
            raise ValueError(f"Unrecognized file type in {file_path}")

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