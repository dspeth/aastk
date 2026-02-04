#!/usr/bin/env python3
from pathlib import Path
from typing import Optional
import shutil
import logging
import zlib
import random
import subprocess

logger = logging.getLogger(__name__)


def bin_mid(bin_series):
	return bin_series.apply(lambda b: (b.left + b.right) / 2)

def check_dependency_availability(command: str):
	if shutil.which(command) is None:
		logger.error(f"Command not found in PATH: {command}")
		raise FileNotFoundError(f"Command not found in PATH: {command}")

def compress_sequence(sequence: str) -> bytes:
	return zlib.compress(sequence.encode('utf-8'), level=9)

def count_fasta_sequences(fasta_path: str) -> int:
	with open(fasta_path, 'r') as f:
		return sum(1 for line in f if line.startswith('>'))

def decompress_sequence(compressed: bytes) -> str:
	return zlib.decompress(compressed).decode('utf-8')

def determine_dataset_name(file: str, splitter: str, part: int, suffix: str = None):
	basename = Path(file).name
	dataset = basename.split(splitter)[part]

	if suffix:
		dataset = dataset.removesuffix(suffix)

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
		if final_path.is_dir():
			pass
		elif force:
			logger.warning(f"Removing existing file: {final_path}")
			final_path.unlink()
		else:
			raise FileExistsError(f"Path already exists: {final_path}")

	if not final_path.exists():
		final_path.parent.mkdir(parents=True, exist_ok=True)
	elif final_path.is_file():
		final_path.parent.mkdir(parents=True, exist_ok=True)

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

def fasta_subsample(fasta: str,
					output_dir: str,
					subset_size: int,
					force: bool = False):
	"""
    Randomly selects subset of pre-defined size from input fasta.

    Args:
        fasta (str): Input fasta file (in most cases homologous dataset as output by PASR)
        output_dir (str): Output directory
        subset_size (int): Number of sequences randomly selected from input fasta
        force (bool): If set, existing files will be overwritten

    Returns:
        output_path: Path to output subset fasta file
    """
	# ===============================
	# Output file path setup
	# ===============================
	dataset = determine_dataset_name(fasta, '.', 0)
	output_file = f"{dataset}_subsample.faa"
	output_path = ensure_path(output_dir, output_file, force=force)

	# ===============================
	# Parse FASTA to dict
	# ===============================
	sequences = read_fasta_to_dict(fasta)
	total_sequences = len(sequences)
	logger.info(f"Total sequences in FASTA: {total_sequences}")

	# ===========================================================================
	# Check if subset size is larger than number of sequences in source FASTA
	# ===========================================================================
	if subset_size > total_sequences:
		logger.warning(
			f"Requested subset size ({subset_size}) is larger than total sequences ({total_sequences}). Using all sequences.")
		subset_size = total_sequences

	# ======================================================
	# Random sampling of seed FASTA for subset creation
	# ======================================================
	subset_keys = random.sample(list(sequences.keys()), subset_size)
	subset_dict = {k: sequences[k] for k in subset_keys}

	# =====================================
	# Write FASTA subset to output file
	# =====================================
	with open(output_path, 'w') as f:
		for header, seq in subset_dict.items():
			f.write(f">{header}\n")
			f.write(f'{seq}\n')

	return output_path

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

def run_diamond_alignment(fasta: str,
						  align_subset: str,
						  subset_size: int,
						  threads: int,
						  blast_columns: list,
						  output: str = None,
						  force: bool = False):
	"""
    Run DIAMOND makedb and blastp to align a full FASTA file to a subset.

    Args:
        fasta (str): Path to the full input FASTA file (query)
        align_subset (str): Path to the subset FASTA file to use as the DIAMOND database (reference)
        subset_size (int): Number of target sequences
        threads (int): Number of threads to use
        output (str): Output directory
        force (bool): If set, existing files will be overwritten

    Returns:
        align_output (str): Path to Blast Tabular Output file for the alignment
    """
	check_dependency_availability('diamond')

	if Path(fasta).is_file():
		pass
	else:
		logger.error("Input seed FASTA not found")
		raise FileNotFoundError(f"FASTA file does not exist: {fasta}")

	logger.info(f"Starting DIAMOND alignment process")
	logger.info(f"Query FASTA: {fasta}")
	logger.info(f"Reference subset: {align_subset}")
	logger.info(f"Using {threads} threads")

	# ===============================
	# Output file path setup
	# ===============================
	dataset = determine_dataset_name(fasta, '.', 0)
	dbname = ensure_path(output, f"{dataset}_subset_dmnd", force=force)
	align_output = ensure_path(output, f"{dataset}_align", force=force)

	# ===============================
	# Subprocess command creation
	# ===============================
	# ==================================
	# Create DIAMOND DB from subset
	# ==================================
	logger.info(f"Creating DIAMOND database: {dbname}")
	run_dmnd_makedb = [
		"diamond", "makedb",
		"--in", align_subset,
		"-d", dbname,
		"-p", str(threads)
	]

	# ===========================================================
	# DIAMOND blastp with seed FASTA as query and subset as DB
	# ===========================================================
	logger.info(f"Running DIAMOND blastp alignment")
	run_dmnd_blastp = [
		"diamond", "blastp",
		"-q", fasta,
		"-d", dbname,
		"-o", align_output,
		"-p", str(threads),
		"-k", str(subset_size),
		"--sensitive",
		"--masking", "0",
		"--comp-based-stats", "0",
		"--outfmt", "6"
	]

	run_dmnd_blastp = run_dmnd_blastp + blast_columns

	# ===============================
	# Run commands created above
	# ===============================
	logger.debug(f"DIAMOND makedb command: {' '.join(run_dmnd_makedb)}")
	try:
		proc = subprocess.Popen(
			run_dmnd_makedb,
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

	db_file = Path(f"{dbname}.dmnd")
	if not db_file.exists():
		logger.error(f"DIAMOND database file not found at {db_file}")
		raise RuntimeError(f"DIAMOND database was not created at {db_file}")

	logger.info(f"Successfully built DIAMOND database at {dbname}")

	logger.debug(f"DIAMOND blastp command: {' '.join(run_dmnd_blastp)}")
	try:
		proc = subprocess.Popen(
			run_dmnd_blastp,
			stdout=subprocess.DEVNULL,
			stderr=subprocess.PIPE,
			text=True
		)
		_, stderr = proc.communicate()

		if proc.returncode != 0:
			logger.error(f"DIAMOND blastp failed with return code {proc.returncode}")
			if stderr:
				logger.error(f"STDERR: {stderr}")
			raise RuntimeError(f"DIAMOND blastp search failed with return code {proc.returncode}")

		if stderr:
			logger.log(99, stderr)
	except FileNotFoundError as e:
		logger.error("DIAMOND executable not found")
		raise RuntimeError("DIAMOND executable not found. Is it installed and in PATH?") from e
	except Exception as e:
		if not isinstance(e, RuntimeError):
			logger.error(f"Unexpected error in DIAMOND blastp search: {e}")
			raise RuntimeError(f"DIAMOND blastp search failed: {e}") from e
		raise

	if not Path(align_output).exists():
		logger.error(f"DIAMOND output file not found at {align_output}")
		raise RuntimeError(f"DIAMOND search did not produce output at {align_output}")

	logger.info(f"Successfully completed DIAMOND search. Results at {align_output}")

	return align_output

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
			yield header, sequence


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
					yield header, sequence

