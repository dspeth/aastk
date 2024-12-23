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
