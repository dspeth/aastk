from .util import *

import sys
import random
import numpy as np
import subprocess
import pandas as pd
import openTSNE
from sklearn.cluster import DBSCAN
import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.DEBUG,  # or logging.INFO if you want less verbosity
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    stream=sys.stdout
)
logging.getLogger('matplotlib').setLevel(logging.WARNING)
logging.getLogger("PIL").setLevel(logging.INFO)

def fasta_subsample(fasta: str, output_path: str, subset_size: int):
    sequences = read_fasta_to_dict(fasta)
    subset_keys = random.sample(list(sequences.keys()), subset_size)
    subset_dict = {k: sequences[k] for k in subset_keys}
    with open(output_path, 'w') as f:
        for header, seq in subset_dict.items():
            f.write(f">{header}\n")
            f.write(f'{seq}\n')

    return output_path

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
    early_tsv = early_df.to_csv(f"{basename}_tsne_early_clust.tsv", sep="\t", index=False)

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

    final_tsv = final_df.to_csv(f"{basename}_tsne_embed.tsv", sep="\t", index=False)


    return early_tsv, final_tsv



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
              metadata_genome: str = None,
              subset: bool = False):
    if subset and not subset_size:
        logger.error("Subset size required if subset is True")
        exit()

    if subset:
        subset_fasta = fasta_subsample(fasta, seed_fasta, subset_size)
    else:
        subset_fasta = seed_fasta
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

def asm_plot(tsv_file: str,
             output: str):
    df = pd.read_csv(tsv_file, sep='\t')

    plt.figure(figsize=(12, 8))

    noise_mask = df['cluster'] == -1
    clustered_mask = df['cluster'] != -1

    if noise_mask.any():
        plt.scatter(df.loc[noise_mask, 'tsne1'],
                    df.loc[noise_mask, 'tsne2'],
                    c='lightgray',
                    alpha=0.6,
                    s=20,
                    label='Noise')

    if clustered_mask.any():
        unique_clusters = df.loc[clustered_mask, 'cluster'].unique()
        colors = plt.cm.Set3(np.linspace(0, 1, len(unique_clusters)))

        for i, cluster in enumerate(unique_clusters):
            cluster_mask = df['cluster'] == cluster
            plt.scatter(df.loc[cluster_mask, 'tsne1'],
                        df.loc[cluster_mask, 'tsne2'],
                        c=[colors[i]], s=30, alpha=0.7,
                        label=f'Cluster {cluster}')

    plt.xlabel('t-SNE 1')
    plt.ylabel('t-SNE 2')
    plt.title('t-SNE Embedding with DBSCAN Clusters')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()

    if output:
        plt.savefig(f"{output}_tsne_clusters.png", dpi=300, bbox_inches='tight')
    plt.show()



