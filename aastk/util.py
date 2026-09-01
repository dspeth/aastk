#!/usr/bin/env python3
from pathlib import Path
from typing import Optional
import shutil
import logging
import zlib
import random
import subprocess
import sqlite3
import json
import yaml
from tqdm import tqdm
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import pandas as pd

from aastk.db.schema import ANNOTATION_COLUMNS


logger = logging.getLogger(__name__)

FILTER_BLAST_OUTPUT_COLUMNS = ['qseqid', 'sseqid', 'nident', 'length', 'qlen']

CUGO_CONTEXT_BASE_COLUMNS = ['target_id', 'position', 'seqID', 'parent_ID',
							 'aa_length', 'strand', 'CUGO_number', 'no_TMH']

BSR_TSV_COLUMNS = ["qseqid", "sseqid", "pident", "qlen", "score", "max_score", "BSR"]


def bin_mid(bin_series):
	return bin_series.apply(lambda b: (b.left + b.right) / 2)

def check_dependency_availability(command: str):
	if shutil.which(command) is None:
		logger.error(f"Command not found in PATH: {command}")
		raise FileNotFoundError(f"Command not found in PATH: {command}")

def compress_sequence(sequence: str) -> bytes:
	return zlib.compress(sequence.encode('utf-8'), level=9)

def count_fasta_sequences(fasta_path: str) -> int:
	if not is_fasta(fasta_path):
		raise ValueError(f"Not a valid FASTA file: {fasta_path}")
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

def ensure_path(path: Optional[str] = None,
				target: Optional[str] = None,
				suffix: Optional[str] = None,
				force: bool = False):
	path = Path(path) if path else Path('.')
	final_path = path / target if target else path

	# add suffix to final_path for path check and optional removal
	if suffix:
		final_path = Path(str(final_path) + f'{suffix}')

	if final_path.exists():
		if final_path.is_dir():
			pass
		elif force:
			logger.warning(f"Removing existing file: {final_path}")
			final_path.unlink()
		else:
			raise FileExistsError(f"Path already exists: {final_path}")

	# remove suffix again
	if suffix:
		final_path = Path(str(final_path).replace(suffix, ''))

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

def filter(fasta: str,
           db_path: str,
           output: str,
           threads: int,
           sql: bool = False,
           svg: bool = False,
           force: bool = False):
    if sql and not db_path:
        raise ValueError('SQL mode requires db_path')

    prefix = determine_dataset_name(fasta, '.', 0)
    output_path = ensure_path(output, f'{prefix}_filtered.faa', force=force)

    if svg:
        plot_path = ensure_path(output, f'{prefix}_filtered.svg', force=force)
    else:
        plot_path = ensure_path(output, f'{prefix}_filtered.png', force=force)

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3)

    subset = fasta_subsample(fasta, output, 100, force=force)

    align_output = run_diamond_alignment(fasta, subset, None, threads, FILTER_BLAST_OUTPUT_COLUMNS, output, force)

    alignment_df = pd.read_csv(align_output, sep='\t', names=FILTER_BLAST_OUTPUT_COLUMNS)

    alignment_df['unaligned_length'] = alignment_df['qlen'] - alignment_df['length']

    means = alignment_df.groupby('qseqid').mean(numeric_only=True)
    means.rename(columns={'nident': 'mean100_nident', 'length': 'mean100_length', 'qlen': 'mean100_qlen', 'unaligned_length': 'mean100_unaligned_length'}, inplace=True)

    mean_avg_length = means.loc[:, 'mean100_length'].mean()

    # first pass filter histogram
    binwidth = 10
    ax1.hist(means['mean100_length'], bins=range(round(min(means['mean100_length'])), round(max(means['mean100_length'])), binwidth))
    ax1.axvline(x=(mean_avg_length + 150), color='black', linestyle='dashed', linewidth=1)
    ax1.axvline(x=(mean_avg_length - 150), color='black', linestyle='dashed', linewidth=1)
    ax1.set_xlabel('avg. alignment length')
    ax1.set_title(label='First pass:\nmean avg. align. length\n+/- 150', fontdict={'fontsize': 10})

    count = 0
    for qseqid in means.index:
        avg_length = means.loc[qseqid, 'mean100_length']
        avg_length_deviation = avg_length - mean_avg_length
        if abs(avg_length_deviation) >= 150:
            means.drop(qseqid, inplace=True)
            count += 1



    remaining = len(means.index)
    logger.info(f"First pass: dropped {count} sequences. Remaining sequences: {remaining}")

    updated_mean_avg_length = means.loc[:, 'mean100_length'].mean()
    updated_std_avg_length = means.loc[:, 'mean100_length'].std()
    lower_bound = updated_mean_avg_length - 3 * updated_std_avg_length
    upper_bound = updated_mean_avg_length + 3 * updated_std_avg_length

    # second pass filter histogram
    ax2.hist(means['mean100_length'], bins=range(round(min(means['mean100_length'])), round(max(means['mean100_length'])), binwidth))
    ax2.axvline(x=lower_bound, color='black', linestyle='dashed', linewidth=1)
    ax2.axvline(x=upper_bound, color='black', linestyle='dashed', linewidth=1)
    ax2.set_xlabel('avg. alignment length')
    ax2.set_title('Second pass:\nmean avg. align. length\n+/- 3 SD', fontdict={'fontsize': 10})

    count = 0
    for qseqid in means.index:
        avg_length = means.loc[qseqid, 'mean100_length']

        if avg_length < lower_bound or avg_length > upper_bound:
            means.drop(qseqid, inplace=True)
            count += 1

    remaining = len(means.index)
    logger.info(f"Second pass: dropped {count} sequences. Remaining sequences: {remaining}")

    penultimate_mean_avg_length = means.loc[:, 'mean100_length'].mean()
    boundary = 0.5 * penultimate_mean_avg_length

    ax3.hist(means['mean100_unaligned_length'],
             bins=range(round(min(means['mean100_unaligned_length'])), round(max(means['mean100_unaligned_length'])), binwidth))
    ax3.axvline(x=boundary, color='black', linestyle='dashed', linewidth=1)
    ax3.set_xlabel("avg. unaligned length")
    ax3.set_title('Third pass:\nmean unaligned length >\n0.5 * mean align. length', fontdict={'fontsize': 10})

    count = 0
    for qseqid in means.index:
        mean_unaligned_length = means.loc[qseqid, 'mean100_unaligned_length']
        boundary = 0.5 * penultimate_mean_avg_length

        if abs(mean_unaligned_length) > boundary:
            means.drop(qseqid, inplace=True)
            count += 1

    remaining = len(means.index)
    logger.info(f"Third pass: dropped {count} sequences. Remaining sequences: {remaining}")

    seq_ids = means.index.dropna().unique().tolist()
    seq_ids = [str(seq_id) for seq_id in seq_ids]

    if sql:
        _ = retrieve_sequences_from_db(seq_ids, output_path, db_path)
    else:
        sequences_written = 0
        with open(output_path, 'w') as f:
            for header, sequence in write_fa_matches(fasta, seq_ids):
                f.write(f"{header}\n{sequence}\n")
                sequences_written += 1
        logger.info(f"Retrieved {sequences_written} sequences to {output_path}")

    plt.subplots_adjust(wspace=1.2)
    #fig.suptitle('Distribution of sequences during filtering with cutoffs')
    plt.savefig(plot_path, dpi=300)

    return output_path

def is_blast_tab(file_path: str, n_columns: int, numeric_columns: list = None) -> bool:
	with open(file_path, 'r') as file:
		first_line = file.readline()
	if not first_line or first_line.startswith('>') or first_line.startswith('@'):
		return False
	fields = first_line.rstrip('\n\r').split('\t')
	if len(fields) != n_columns or any(field == '' for field in fields):
		return False
	if numeric_columns:
		for idx in numeric_columns:
			try:
				float(fields[idx])
			except ValueError:
				return False
	return True

def is_bsr_tsv(file_path: str) -> bool:
	with open(file_path, 'r') as f:
		header = f.readline().rstrip('\n\r').split('\t')
		first_data_line = f.readline()
	if header != BSR_TSV_COLUMNS:
		return False
	if first_data_line:
		parts = first_data_line.rstrip('\n\r').split('\t')
		if len(parts) != len(BSR_TSV_COLUMNS):
			return False
		try:
			float(parts[2])
			int(parts[3])
			float(parts[4])
			float(parts[5])
			float(parts[6])
		except ValueError:
			return False
	return True

def is_casm_embedding_tsv(file_path: str) -> bool:
	with open(file_path, 'r') as f:
		header_line = f.readline().rstrip('\n\r')
		first_data_line = f.readline()
	header = header_line.split('\t')
	if 'seqID' not in header or 'cluster' not in header:
		return False
	if first_data_line:
		parts = first_data_line.rstrip('\n\r').split('\t')
		if len(parts) != len(header):
			return False
		try:
			int(parts[header.index('cluster')])
		except ValueError:
			return False
	return True

def is_casm_matrix_metadata_json(file_path: str) -> bool:
	try:
		with open(file_path, 'r') as f:
			data = json.load(f)
	except (json.JSONDecodeError, OSError):
		return False
	return (isinstance(data, dict)
			and isinstance(data.get('queries'), list)
			and isinstance(data.get('targets'), list))

def is_column_info_json(file_path: str) -> bool:
	try:
		with open(file_path, 'r') as f:
			data = json.load(f)
	except (json.JSONDecodeError, OSError):
		return False
	required = ('score', 'qlen', 'pident')
	return isinstance(data, dict) and all(isinstance(data.get(k), int) for k in required)

def is_cugo_context_tsv(file_path: str) -> bool:
	with open(file_path, 'r') as f:
		header_line = f.readline().rstrip('\n\r')
		first_data_line = f.readline()
	header = header_line.split('\t')
	if not all(col in header for col in CUGO_CONTEXT_BASE_COLUMNS):
		return False
	if not any(col in header for col in ANNOTATION_COLUMNS):
		return False
	if first_data_line:
		parts = first_data_line.rstrip('\n\r').split('\t')
		if len(parts) != len(header):
			return False
		row = dict(zip(header, parts))
		try:
			int(row['position'])
		except ValueError:
			return False
	return True

def is_fasta(file_path) -> bool:
	with open(file_path, 'r') as file:
		return file.read(1) == '>'

def is_fastq(file_path) -> bool:
	with open(file_path, 'r') as file:
		return file.read(1) == '@'

def is_max_score_tsv(file_path: str) -> bool:
	with open(file_path, 'r') as f:
		header = f.readline().rstrip('\n\r')
		first_data_line = f.readline()
	if header != "Protein_id\tmax_score":
		return False
	if first_data_line:
		parts = first_data_line.rstrip('\n\r').split('\t')
		if len(parts) != 2:
			return False
		try:
			float(parts[1])
		except ValueError:
			return False
	return True

def is_pasr_threshold_yaml(file_path: str) -> bool:
	try:
		with open(file_path) as f:
			data = yaml.safe_load(f)
	except (yaml.YAMLError, OSError):
		return False
	return isinstance(data, dict) and 'max_score_min' in data and 'max_score_max' in data

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
	if not is_fasta(fasta):
		raise ValueError(f"Not a valid FASTA file: {fasta}")

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

def retrieve_sequences_from_db(seq_ids: list,
							   output_path: str,
							   db_path: str):
	if not seq_ids:
		logger.warning(f"No sequences found for position {position}")
		return output_path

	# connect to database
	conn = sqlite3.connect(db_path)

	# query database for sequences
	placeholders = ','.join('?' * len(seq_ids))
	cursor = conn.execute(f"""
			SELECT seqID, protein_seq
			FROM protein_data
			WHERE seqID IN ({placeholders}) AND protein_seq IS NOT NULL
		""", seq_ids)

	# write to file
	with open(output_path, 'w') as file:
		count = 0
		for seqid, compressed_seq in tqdm(cursor, total=len(seq_ids), desc="Retrieving sequences"):
			sequence = decompress_sequence(compressed_seq)
			file.write(f">{seqid}\n{sequence}\n")
			count += 1

	conn.close()
	logger.info(f"Retrieved {count} sequences to {output_path}")

	return output_path

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

	if not is_fasta(fasta):
		raise ValueError(f"Not a valid FASTA file: {fasta}")

	if not is_fasta(align_subset):
		raise ValueError(f"Not a valid FASTA file: {align_subset}")

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
	if not is_fasta(seq_file):
		raise ValueError(f"Not a valid FASTA file: {seq_file}")

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
	if not is_fastq(seq_file):
		raise ValueError(f"Not a valid FASTQ file: {seq_file}")

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
