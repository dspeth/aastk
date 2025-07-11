import json

from .util import *

import sys
import random
import numpy as np
import subprocess
import pandas as pd
import openTSNE
from sklearn.cluster import DBSCAN
import matplotlib.pyplot as plt
import yaml
import re

logger = logging.getLogger(__name__)
logging.basicConfig(
	level=logging.DEBUG,  # or logging.INFO if you want less verbosity
	format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
	stream=sys.stdout
)
logging.getLogger('matplotlib').setLevel(logging.WARNING)
logging.getLogger("PIL").setLevel(logging.INFO)


def fasta_subsample(fasta: str,
                    output_dir: str,
                    subset_size: int,
                    force: bool = False):
	logger.info(f"Starting FASTA subsampling from {fasta}")
	dataset = determine_dataset_name(fasta, '.', 0)
	output_file = f"{dataset}_subsample.faa"
	output_path = ensure_path(output_dir, output_file, force=force)

	logger.debug(f"Reading FASTA file: {fasta}")
	sequences = read_fasta_to_dict(fasta)
	total_sequences = len(sequences)
	logger.info(f"Total sequences in FASTA: {total_sequences}")

	if subset_size > total_sequences:
		logger.warning(
			f"Requested subset size ({subset_size}) is larger than total sequences ({total_sequences}). Using all sequences.")
		subset_size = total_sequences

	logger.info(f"Randomly sampling {subset_size} sequences from {total_sequences}")
	subset_keys = random.sample(list(sequences.keys()), subset_size)
	subset_dict = {k: sequences[k] for k in subset_keys}

	logger.debug(f"Writing subset to: {output_path}")
	with open(output_path, 'w') as f:
		for header, seq in subset_dict.items():
			f.write(f">{header}\n")
			f.write(f'{seq}\n')

	logger.info(f"Successfully created subset FASTA with {len(subset_dict)} sequences")
	return output_path


def run_diamond_alignment(fasta: str,
                          align_subset: str,
                          subset_size: int,
                          threads: int):
	"""
    Run DIAMOND makedb and blastp to align a full FASTA file to a subset.

    Args:
        fasta (str): Path to the full input FASTA file (query).
        align_subset (str): Path to the subset FASTA file to use as the DIAMOND database (subject).
        subset_size (int): Number of target sequences (used for -k).
        threads (int): Number of threads to use.
    """
	logger.info(f"Starting DIAMOND alignment process")
	logger.info(f"Query FASTA: {fasta}")
	logger.info(f"Reference subset: {align_subset}")
	logger.info(f"Using {threads} threads")

	dataset = determine_dataset_name(fasta, '.', 0)
	dbname = f"{dataset}_subset_dmnd"
	align_output = f"{dataset}_align"

	logger.info(f"Creating DIAMOND database: {dbname}")
	run_dmnd_makedb = [
		"diamond", "makedb",
		"--in", align_subset,
		"-d", dbname,
		"-p", str(threads)
	]

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
		"--outfmt", "6", "qseqid", "sseqid", "score",
		"--comp-based-stats", "0"
	]

	logger.debug(f"DIAMOND makedb command: {' '.join(run_dmnd_makedb)}")
	subprocess.run(run_dmnd_makedb)
	logger.info("DIAMOND database creation completed")

	logger.debug(f"DIAMOND blastp command: {' '.join(run_dmnd_blastp)}")
	subprocess.run(run_dmnd_blastp)
	logger.info(f"DIAMOND alignment completed. Output saved to: {align_output}")

	return align_output

def build_alignment_matrix_split(align_file: str,
                                 output: str = None,
                                 force: bool = False):
	logger.info(f"Building alignment matrix from: {align_file} (split parsing)")

	queries_set = set()
	targets_set = set()

	with open(align_file, 'r') as f:
		for line_num, line in enumerate(f, 1):
			line = line.rstrip('\n\r')
			if not line:
				continue

			try:
				parts = line.split("\t")
				if len(parts) != 3:
					logger.warning(f"Line {line_num}: Expected 3 fields, got {len(parts)}. Skipping.")
					continue

				query, target, score_str = parts
				queries_set.add(query)
				targets_set.add(target)
			except Exception as e:
				logger.warning(f"Line {line_num}: Error parsing line. Skipping. ({e})")
				continue

	queries = sorted(queries_set)
	targets = sorted(targets_set)
	query_to_idx = {q: i for i, q in enumerate(queries)}
	target_to_idx = {t: i for i, t in enumerate(targets)}

	del queries_set, targets_set

	matrix = np.zeros((len(queries), len(targets)), dtype=np.float32)

	with open(align_file, 'r') as f:
		for line_num, line in enumerate(f, 1):
			line = line.rstrip('\n\r')
			if not line:
				continue

			try:
				parts = line.split("\t")
				if len(parts) != 3:
					continue

				query, target, score_str = parts
				score = float(score_str)
				i = query_to_idx[query]
				j = target_to_idx[target]
				matrix[i, j] = score
			except ValueError as e:
				logger.warning(f"Line {line_num}: Invalid score '{score_str}'. Skipping. ({e})")
				continue
			except KeyError as e:
				logger.warning(f"Line {line_num}: Key not found in mapping. Skipping. ({e})")
				continue
			except Exception as e:
				logger.warning(f"Line {line_num}: Unexpected error. Skipping. ({e})")
				continue

	non_zero_elements = np.count_nonzero(matrix)
	matrix_elements = len(queries) * len(targets)
	sparsity = (matrix_elements - non_zero_elements) / matrix_elements * 100

	logger.info(f"Matrix construction completed")
	logger.info(f"Non-zero elements: {non_zero_elements:,} ({100 - sparsity:.2f}% filled)")
	logger.info(f"Matrix sparsity: {sparsity:.2f}%")
	logger.info(f"Score range: {matrix.min():.2f} - {matrix.max():.2f}")

	dataset_name = determine_dataset_name(align_file, '.', 0)

	if output:
		matrix_file = ensure_path(output, f"{dataset_name}_matrix.npy", force=force)
		metadata_file = ensure_path(output, f"{dataset_name}_matrix_metadata.json", force=force)
	else:
		matrix_file = ensure_path(target=f"{dataset_name}_matrix.npy", force=force)
		metadata_file = ensure_path(target=f"{dataset_name}_matrix_metadata.json", force=force)

	np.save(matrix_file, matrix)
	logger.info(f"Matrix saved to: {matrix_file}")

	metadata = {
		"queries": queries,
		"targets": targets,
		"matrix_shape": matrix.shape,
		"non_zero_elements": int(non_zero_elements),
		"sparsity": float(sparsity),
		"score_range": [float(matrix.min()), float(matrix.max())]
	}

	with open(metadata_file, 'w') as file:
		json.dump(metadata, file, indent=2)
	logger.info(f"Matrix metadata saved to: {metadata_file}")

	return matrix, queries, targets

def load_alignment_matrix_from_file(matrix_path: str,
                                    metadata_path: str):
	matrix = np.load(matrix_path)

	with open(metadata_path, "r") as file:
		metadata = json.load(file)

	queries = metadata["queries"]
	targets = metadata["targets"]

	logger.info(f"Loaded matrix shape: {matrix.shape}")
	logger.info(f"Loaded {len(queries)} queries and {len(targets)} targets")

	return matrix, queries, targets

def create_embedding_dataframe(embedding: np.ndarray,
                               queries: list,
                               clusters: np.ndarray,
                               col_names: list,
                               metadata_protein: str = None,
                               metadata_genome: str = None
                               ):
	logger.debug(f"Creating embedding dataframe with {len(queries)} proteins")
	df = pd.DataFrame(embedding, columns=col_names)
	df['prot_ID'] = queries
	df['cluster'] = clusters

	logger.debug("Extracting locus and genome information from protein IDs")
	# extract locus and genome info from protein IDs
	df['locus_nr'] = df['prot_ID'].astype(str).str.rsplit("_", n=1).str[1]
	df['genome_ID'] = df['prot_ID'].astype(str).str.rsplit("_", n=1).str[0]

	# reorder columns
	base_cols = ['prot_ID', 'locus_nr', 'genome_ID'] + col_names + ['cluster']
	df = df[base_cols]

	# add metadata if provided
	if metadata_protein and metadata_protein.is_file():
		logger.info(f"Adding protein metadata from: {metadata_protein}")
		protein_meta = pd.read_csv(metadata_protein, sep="\t")
		df = pd.merge(df, protein_meta, on="prot_ID", how="left")

	if metadata_genome and metadata_genome.is_file():
		logger.info(f"Adding genome metadata from: {metadata_genome}")
		genome_meta = pd.read_csv(metadata_genome, sep="\t")
		df = pd.merge(df, genome_meta, on="genome_ID", how="left")

	logger.debug(f"Final dataframe shape: {df.shape}")
	return df


def tsne_embedding(matrix: np.ndarray,
                   queries: list,
                   basename: str,
                   perplexity: int,
                   iterations: int,
                   exaggeration: int,
                   threads: int,
                   metadata_protein: str = None,
                   metadata_genome: str = None,
                   force: bool = False
                   ):
	logger.info("Starting t-SNE embedding process")
	logger.info(f"Matrix shape: {matrix.shape}")
	logger.info(f"Perplexity: {perplexity}, Iterations: {iterations}, Exaggeration: {exaggeration}")
	logger.info(f"Using {threads} threads")

	# t-SNE embedding
	logger.info("Computing multiscale affinities")
	affinities = openTSNE.affinity.Multiscale(matrix,
	                                          perplexities=[20, perplexity],
	                                          metric="cosine",
	                                          n_jobs=threads)

	logger.info("Initializing t-SNE with PCA")
	pca_init = openTSNE.initialization.pca(matrix)
	tsne_embed = openTSNE.TSNEEmbedding(pca_init,
	                                    affinities,
	                                    n_jobs=threads,
	                                    dof=1,
	                                    learning_rate="auto")

	# early embed with exaggeration
	logger.info(f"Starting early optimization with exaggeration={exaggeration}")
	tsne_embed.optimize(n_iter=iterations,
	                    exaggeration=exaggeration,
	                    momentum=0.5,
	                    inplace=True,
	                    n_jobs=threads)
	logger.info("Early optimization completed")

	# cluster on early embedding
	logger.info("Performing DBSCAN clustering on early embedding")
	clustering = DBSCAN(eps=1, min_samples=10).fit(tsne_embed)
	n_clusters = len(set(clustering.labels_)) - (1 if -1 in clustering.labels_ else 0)
	n_noise = list(clustering.labels_).count(-1)
	logger.info(f"DBSCAN clustering results: {n_clusters} clusters, {n_noise} noise points")

	# save early embedding
	logger.debug("Saving early embedding results")
	early_df = create_embedding_dataframe(tsne_embed,
	                                      queries,
	                                      clustering.labels_,
	                                      ['tsne1', 'tsne2'],
	                                      metadata_protein,
	                                      metadata_genome,
	                                      )
	early_filename = ensure_path(target=f"{basename}_tsne_early_clust.tsv", force=force)
	early_tsv = early_df.to_csv(early_filename, sep="\t", index=False)
	logger.info(f"Early embedding saved to: {early_filename}")

	# final embedding
	logger.info("Starting final optimization")
	tsne_embed.optimize(n_iter=iterations,
	                    momentum=0.8,
	                    inplace=True,
	                    n_jobs=threads)
	logger.info("Final optimization completed")

	# save final embedding
	logger.debug("Saving final embedding results")
	final_df = create_embedding_dataframe(tsne_embed,
	                                      queries,
	                                      clustering.labels_,
	                                      ['tsne1', 'tsne2'],
	                                      metadata_protein,
	                                      metadata_genome,
	                                      )

	final_filename = ensure_path(target=f"{basename}_tsne_embed.tsv", force=force)
	final_tsv = final_df.to_csv(final_filename, sep="\t", index=False)
	logger.info(f"Final embedding saved to: {final_filename}")

	return early_tsv, final_tsv, early_filename,final_filename

def casm_plot(tsv_file: str,
              output: str,
              ):
	logger.info(f"Creating t-SNE plot from: {tsv_file}")
	df = pd.read_csv(tsv_file, sep='\t')
	logger.info(f"Loaded {len(df)} data points for plotting")

	plt.figure(figsize=(14, 10))

	noise_mask = df['cluster'] == -1
	clustered_mask = df['cluster'] != -1

	n_noise = noise_mask.sum()
	n_clustered = clustered_mask.sum()
	logger.info(f"Plotting {n_clustered} clustered points and {n_noise} noise points")

	if noise_mask.any():
		plt.scatter(df.loc[noise_mask, 'tsne1'],
		            df.loc[noise_mask, 'tsne2'],
		            c='lightgray',
		            alpha=0.4,
		            s=12,
		            edgecolors='none',
		            label='Noise')

	if clustered_mask.any():
		unique_clusters = df.loc[clustered_mask, 'cluster'].unique()
		colors = plt.cm.tab20(np.linspace(0, 1, len(unique_clusters)))
		if len(unique_clusters) > 20:
			colors = plt.cm.hsv(np.linspace(0, 1, len(unique_clusters)))

		logger.info(f"Plotting {len(unique_clusters)} unique clusters")

		for i, cluster in enumerate(unique_clusters):
			cluster_mask = df['cluster'] == cluster
			cluster_size = cluster_mask.sum()
			logger.debug(f"Cluster {cluster}: {cluster_size} points")
			plt.scatter(df.loc[cluster_mask, 'tsne1'],
			            df.loc[cluster_mask, 'tsne2'],
			            c=[colors[i]],
			            s=15,
			            alpha=0.7,
			            edgecolors='white',
			            linewidths=0.3,
			            label=f'Cluster {cluster}')

	plt.xlabel('t-SNE 1', fontsize=12)
	plt.ylabel('t-SNE 2', fontsize=12)
	plt.title('t-SNE Embedding with DBSCAN Clusters', fontsize=14)

	if len(df.loc[clustered_mask, 'cluster'].unique()) <= 15:
		plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)

	plt.tight_layout()

	if output:
		suffix = "_tsne_cog.png" if 'cog' in locals() else "_tsne_clusters.png"
		output_file = f"{output}{suffix}"
		plt.savefig(output_file, dpi=300, bbox_inches='tight')
		logger.info(f"Plot saved to: {output_file}")
	plt.show()

def casm(fasta: str,
         output: str,
         subset: str,
         subset_size: int,
         threads: int = 1,
         perplexity: int = 50,
         iterations: int = 500,
         exaggeration: int = 6,
         metadata_protein: str = None,
         metadata_genome: str = None,
         matrix_path: str = None,
         matrix_metadata_path: str = None,
         early_clust_path: str = None,
         full_clust_path: str = None,
         force: bool = False,
         matrix: bool = False,
         embed: bool = False,
         plot: bool = False,
         all: bool = False
         ):
	if matrix and not fasta:
		logger.error("No input fasta provided. Alignment matrix can not be created.")
		exit()

	if embed:
		if not matrix_path:
			logger.error("No input matrix file provided. No TSNE embedding possible")
			exit()
		if not matrix_metadata_path:
			logger.error("No input matrix metadata file provided. No TSNE embedding possible")
			exit()

	if plot:
		if not early_clust_path:
			logger.error("No input early clustering TSV provided. Plotting not possible.")
			exit()
		if not full_clust_path:
			logger.error("No input full clustering TSV provided. Plotting not possible.")
			exit()

	if all and not fasta:
		logger.error("No input fasta provided. Full casm workflow not possible.")
		exit()

	if all:
		logger.info("=== Starting CASM Analysis ===")
		logger.info(f"Input FASTA: {fasta}")
		logger.info(f"Output basename: {output}")
		logger.info(f"Subset size: {subset_size}")
		logger.info(f"Threads: {threads}")

		if subset:
			logger.info(f"Using provided subset: {subset}")
			subset_fasta = subset
		else:
			logger.info("Creating random subset from input FASTA")
			subset_fasta = fasta_subsample(fasta, output, subset_size, force=force)

		logger.info("=== Phase 1: DIAMOND Alignment ===")
		align_output = run_diamond_alignment(fasta, subset_fasta, subset_size, threads)

		logger.info("=== Phase 2: Matrix Construction ===")
		matrix, queries, targets = build_alignment_matrix_split(align_output, output, force)

		logger.info("=== Phase 3: t-SNE Embedding ===")
		_, _, early_filename, final_filename = tsne_embedding(matrix=matrix,
		                                                      queries=queries,
		                                                      basename=output,
		                                                      perplexity=perplexity,
		                                                      iterations=iterations,
		                                                      exaggeration=exaggeration,
		                                                      threads=threads,
		                                                      metadata_protein=metadata_protein,
		                                                      metadata_genome=metadata_genome,
		                                                      force=force
		                                                      )

		logger.info("=== Phase 4: Generating plots for early and full clustering ===")
		casm_plot(early_filename, output=output)
		casm_plot(final_filename, output=output)

		logger.info("=== CASM Analysis Completed ===")

	elif matrix:
		logger.info("=== Starting Matrix Construction ===")
		logger.info(f"Input FASTA: {fasta}")
		logger.info(f"Output basename: {output}")
		logger.info(f"Subset size: {subset_size}")
		logger.info(f"Threads: {threads}")

		if subset:
			logger.info(f"Using provided subset: {subset}")
			subset_fasta = subset
		else:
			logger.info("Creating random subset from input FASTA")
			subset_fasta = fasta_subsample(fasta, output, subset_size, force=force)

		logger.info("=== Phase 1: DIAMOND Alignment ===")
		align_output = run_diamond_alignment(fasta, subset_fasta, subset_size, threads)

		logger.info("=== Phase 2: Matrix Construction ===")
		build_alignment_matrix_split(align_output, output, force)

	elif embed:
		logger.info("=== t-SNE Embedding ===")
		matrix, queries, targets = load_alignment_matrix_from_file(matrix_path, matrix_metadata_path)
		tsne_embedding(matrix=matrix,
                       queries=queries,
                       basename=output,
                       perplexity=perplexity,
                       iterations=iterations,
                       exaggeration=exaggeration,
                       threads=threads,
                       metadata_protein=metadata_protein,
                       metadata_genome=metadata_genome,
		               force=force
                       )

	elif plot:
		logger.info("=== Generating plots for early and full clustering ===")
		casm_plot(early_clust_path, output=output)
		casm_plot(full_clust_path, output=output)

