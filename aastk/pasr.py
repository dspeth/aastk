#!/usr/bin/env python3

from .util import extract_unique_keys, determine_file_type

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
            with open(seq_file, "r") as sf:
                for header, sequence in write_fa_matches(sf, matching_ids):
                    out.write(f">{header}\n{sequence}\n")
        elif file_type == "fastq":
            with open(seq_file, "r") as sf:
                for header, sequence in write_fq_matches(sf, matching_ids):
                    out.write(f">{header}\n{sequence}\n")



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
