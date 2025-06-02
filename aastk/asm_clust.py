from .util import *

import sys
from pathlib import Path
import random
import numpy as np
import subprocess
import pandas as pd
import openTSNE
from sklearn.cluster import DBSCAN

def fasta_subsample(fasta: str, output_path: str, subset_size: int):
    sequences = read_fasta_to_dict(fasta)
    subset_keys = random.sample(list(sequences.keys()), subset_size)
    subset_dict = {k: sequences[k] for k in subset_keys}
    with open(output_path, 'w') as f:
        for header, seq in subset_dict.items():
            f.write(f">{header}\n")
            f.write(f'{seq}\n')

    return output_path

import subprocess
from pathlib import Path

def run_diamond_alignment(fasta: str, align_subset: str, dataset: str, subset_size: int, threads: int):
    """
    Run DIAMOND makedb and blastp to align a full FASTA file to a subset.

    Args:
        fasta (str): Path to the full input FASTA file (query).
        align_subset (str): Path to the subset FASTA file to use as the DIAMOND database (subject).
        output_prefix (str): Prefix for output files (e.g., alignment and database).
        subset_size (int): Number of target sequences (used for -k).
        threads (int): Number of threads to use.
    """
    dbname = f"{dataset}_subset_dmnd"
    align_output = f"{dataset}_align"

    run_dmnd_makedb = [
        "diamond", "makedb",
        "--in", align_subset,
        "-d", dbname,
        "-p", str(threads)
    ]

    run_dmnd_blastp = [
        "diamond", "blastp",
        "-q", fasta,
        "-d", dbname,
        "-o", align_output,
        "-p", str(threads),
        "-k", str(subset_size),
        "--sensitive",
        "--masking", "0",
        "--outfmt", "6", "qseqid", "sseqid", "score",
        "--comp-based-stats", "0"
    ]

    subprocess.run(run_dmnd_makedb)
    subprocess.run(run_dmnd_blastp)

    return align_output


def build_alignment_matrix(align_file: str):
    align_data = pd.read_csv(align_file, sep="\t", header=None, names=['query', 'target', 'score'])

    # ordered list of protein IDs as rows
    queries = sorted(align_data['query'].unique())
    # ordered list of reference subset IDs as columns
    targets = sorted(align_data['target'].unique())

    query_to_idx = {q: i for i, q in enumerate(queries)}
    target_to_idx = {t: i for i, t in enumerate(targets)}

    matrix = np.zeros((len(queries), len(targets)), dtype=np.float32)

    for _, row in align_data.iterrows():
        i = query_to_idx[row['query']]
        j = target_to_idx[row['target']]
        matrix[i, j] = row['score']

    return matrix, queries, targets



# subset should be optional, could also be any fasta file as input
def asm_clust(fasta: str,
              output: str,
              subset_size: int,
              seed_fasta: str,
              dataset: str,
              threads: int = 1):
    subset_fasta = fasta_subsample(fasta, seed_fasta, subset_size)
    align_output = run_diamond_alignment(fasta, seed_fasta, dataset, subset_size, threads)
    matrix, query, targets = build_alignment_matrix(align_output)
