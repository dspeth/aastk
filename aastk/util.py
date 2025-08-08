#!/usr/bin/env python3
from pathlib import Path
from typing import Optional
import shutil
import logging

logger = logging.getLogger(__name__)


def bin_mid(bin_series):
	return bin_series.apply(lambda b: (b.left + b.right) / 2)

def check_dependency_availability(command: str):
	if shutil.which(command) is None:
		logger.error(f"Command not found in PATH: {command}")
		raise FileNotFoundError(f"Command not found in PATH: {command}")


def count_fasta_sequences(fasta_path: str) -> int:
	with open(fasta_path, 'r') as f:
		return sum(1 for line in f if line.startswith('>'))

def determine_dataset_name(file: str, splitter:str, part: int):
	dataset = file.split(splitter)[part]
	return dataset

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

def ensure_path(path: Optional[str] = None, target: Optional[str] = None, force: bool = False):
	path = Path(path) if path else Path('.')
	final_path = path / target if target else path

	if final_path.exists():
		if force:
			if final_path.is_dir():
				logger.warning(f"Removing existing directory: {final_path}")
				shutil.rmtree(final_path)
			else:
				logger.warning(f"Removing existing file: {final_path}")
				final_path.unlink()
		else:
			raise FileExistsError(f"Path already exists: {final_path}")

	path.mkdir(parents=True, exist_ok=True)

	return str(final_path)

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

def parse_protein_identifier(id: str):
	"""
	Parse protein identifier to retrieve genome identifier

	Args:
		id (str): protein ID

	Returns:
		genome_id (str): genome ID
	"""
	genome_id = id.rsplit("___", 1)[0]
	return genome_id


def read_fasta_to_dict(fasta: str):
	"""
	Reads a FASTA file and returns a dictionary of headers and sequences.

	Args:
		fasta_path (str): Path to the FASTA file

	Returns:
		dict: Dictionary with headers as keys and sequences as values
	"""
	# Load matched FASTA
	sequences = {}
	current_header = None

	# maybe replace with funct from utils
	with open(fasta) as f:
		for line in f:
			line = line.strip()
			if not line:
				continue
			if line.startswith(">"):
				current_header = line[1:]
				sequences[current_header] = ""
			else:
				if current_header:
					sequences[current_header] += line

	return sequences

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

