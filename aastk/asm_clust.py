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


def create_embedding_dataframe(embedding: np.ndarray,
                               queries: list,
                               clusters: np.ndarray,
                               col_names: list,
                               metadata_protein: str = None,
                               metadata_genome: str = None):
    df = pd.DataFrame(embedding, columns=col_names)
    df['prot_ID'] = queries
    df['cluster'] = clusters

    # extract locus and genome info from protein IDs
    df['locus_nr'] = df['prot_ID'].astype(str).str.rsplit("_", n=1).str[1]
    df['genome_ID'] = df['prot_ID'].astype(str).str.rsplit("_", n=1).str[0]

    # reorder columns
    base_cols = ['prot_ID', 'locus_nr', 'genome_ID'] + col_names + ['cluster']
    df = df[base_cols]

    # add metadata if provided
    if metadata_protein and metadata_protein.is_file():
        protein_meta = pd.read_csv(metadata_protein, sep="\t")
        df = pd.merge(df, protein_meta, on="prot_ID", how="left")

    if metadata_genome and metadata_genome.is_file():
        genome_meta = pd.read_csv(metadata_genome, sep="\t")
        df = pd.merge(df, genome_meta, on="genome_ID", how="left")

    return df

def tsne_embedding(matrix: np.ndarray,
                   queries: list,
                   basename: str,
                   perplexity:int,
                   iterations: int,
                   exaggeration: int,
                   threads: int,
                   metadata_protein: str = None,
                   metadata_genome: str = None):
    # t-SNE embedding
    affinities = openTSNE.affinity.Multiscale(matrix,
                                              perplexities=[20, perplexity],
                                              metric="cosine",
                                              n_jobs=threads)

    pca_init = openTSNE.initialization.pca(matrix)
    tsne_embed = openTSNE.TSNEEmbedding(pca_init,
                                        affinities,
                                        n_jobs=threads,
                                        dof=1,
                                        learning_rate="auto")

    # early embed with exaggeration
    tsne_embed.optimize(n_iter=iterations,
                        exaggeration=exaggeration,
                        momentum=0.5,
                        inplace=True,
                        n_jobs=threads)

    # cluster on early embedding
    clustering = DBSCAN(eps=1, min_samples=10).fit(tsne_embed)

    # save early embedding
    early_df = create_embedding_dataframe(tsne_embed,
                                          queries,
                                          clustering.labels_,
                                          ['tsne1', 'tsne2'],
                                          metadata_protein,
                                          metadata_genome)
    early_df.to_csv(f"{basename}_tsne_early_clust.tsv", sep="\t", index=False)

    # final embedding
    tsne_embed.optimize(n_iter=iterations,
                        momentum=0.8,
                        inplace=True,
                        n_jobs=threads)

    # save final embedding
    final_df = create_embedding_dataframe(tsne_embed,
                                          queries,
                                          clustering.labels_,
                                          ['tsne1', 'tsne2'],
                                          metadata_protein,
                                          metadata_genome)

    final_df.to_csv(f"{basename}_tsne_embed.tsv", sep="\t", index=False)


# subset should be optional, could also be any fasta file as input
def asm_clust(fasta: str,
              output: str,
              subset_size: int,
              seed_fasta: str,
              dataset: str,
              threads: int = 1,
              perplexity: int = 50,
              iterations: int = 500,
              exaggeration: int = 6,
              metadata_protein: str = None,
              metadata_genome: str = None):
    subset_fasta = fasta_subsample(fasta, seed_fasta, subset_size)
    align_output = run_diamond_alignment(fasta, subset_fasta, dataset, subset_size, threads)
    matrix, queries, targets = build_alignment_matrix(align_output)
    tsne_embedding(matrix=matrix,
                   queries=queries,
                   basename=output,
                   perplexity=perplexity,
                   iterations=iterations,
                   exaggeration=exaggeration,
                   threads=threads,
                   metadata_protein=metadata_protein,
                   metadata_genome=metadata_genome)

