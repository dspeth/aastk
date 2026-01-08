from .util import *

import logging
import random
import numpy as np
import subprocess
import pandas as pd
import openTSNE
from sklearn.cluster import DBSCAN
import matplotlib.pyplot as plt
import sqlite3
import json


logger = logging.getLogger(__name__)

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

def run_diamond_alignment(fasta: str,
                          align_subset: str,
                          subset_size: int,
                          threads: int,
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
        "--outfmt", "6", "qseqid", "sseqid", "score",
        "--comp-based-stats", "0"
    ]

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

        if stderr:
            logger.log(99, stderr)

    except subprocess.CalledProcessError as e:
        logger.error(f"Error in building the DIAMOND database: {e}")
        logger.error(f"STDERR: {e.stderr}")
        raise RuntimeError(f"DIAMOND database creation failed: {e}") from e

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

        if stderr:
            logger.log(99, stderr)
    except subprocess.CalledProcessError as e:

        logger.error(f"Error in DIAMOND blastp search: {e}")
        logger.error(f"STDERR: {e.stderr}")
        raise RuntimeError(f"DIAMOND blastp search failed: {e}") from e

    logger.info(f"Successfully completed DIAMOND search. Results at {align_output}")

    return align_output

def build_alignment_matrix_split(align_file: str,
                                 output: str = None,
                                 force: bool = False):
    """
    Builds alignment matrix from input Blast Tabular Output to be used in tSNE embedding

    Args:
        align_file (str): Path to Blast Tabular Output
        output (str): Path to desired output directory
        force (bool): If set, existing files will be overwritten

    Returns:
        matrix: Alignment matrix as dense numpy matrix
        queries (list): List of query protein IDs
        targets (list): List of reference protein IDs
        matrix_file (str): Path to .npy file containing numpy matrix
        metadata_file (str): Path to .json file containing matrix metadata (queries, targets, matrix stats)
    """
    logger.info(f"Building alignment matrix from: {align_file} (split parsing)")

    # ===============================
    # Storage preparation
    # ===============================

    # set up empty sets for storing query and reference protein IDs
    queries_set = set()
    targets_set = set()

    # =======================================================
    # First pass over alignment file: Retrieve protein IDs
    # =======================================================
    with open(align_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.rstrip('\n\r')
            if not line:
                continue

            try:
                # split tab separated line to get query ID, target ID and score
                parts = line.split("\t")
                if len(parts) != 3:
                    logger.warning(f"Line {line_num}: Expected 3 fields, got {len(parts)}. Skipping.")
                    continue
                query, target, score_str = parts

                # add protein IDs to respective sets
                queries_set.add(query)
                targets_set.add(target)
            except Exception as e:
                logger.warning(f"Line {line_num}: Error parsing line. Skipping. ({e})")
                continue

    # sort both sets
    queries = sorted(queries_set)
    targets = sorted(targets_set)

    # determine indices for each ID for matrix construction and save {protID: index}
    query_to_idx = {q: i for i, q in enumerate(queries)}
    target_to_idx = {t: i for i, t in enumerate(targets)}

    # delete ID sets for efficient memory handling
    del queries_set, targets_set

    # set up empty numpy matrix with correct dimensions
    matrix = np.zeros((len(queries), len(targets)), dtype=np.float32)

    # =======================================================
    # Second pass over alignment file: Matrix construction
    # =======================================================
    with open(align_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.rstrip('\n\r')
            if not line:
                continue

            # split tab separated line to get query ID, target ID and score
            try:
                parts = line.split("\t")
                if len(parts) != 3:
                    continue
                query, target, score_str = parts

                score = float(score_str)
                i = query_to_idx[query]
                j = target_to_idx[target]

                # assign alignment score to correct matrix cell
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

    # ===============================
    # Determine matrix statistics
    # ===============================
    non_zero_elements = np.count_nonzero(matrix)
    matrix_elements = len(queries) * len(targets)
    sparsity = (matrix_elements - non_zero_elements) / matrix_elements * 100

    logger.info(f"Matrix construction completed")
    logger.info(f"Non-zero elements: {non_zero_elements:,} ({100 - sparsity:.2f}% filled)")
    logger.info(f"Matrix sparsity: {sparsity:.2f}%")
    logger.info(f"Score range: {matrix.min():.2f} - {matrix.max():.2f}")

    dataset_name = determine_dataset_name(align_file, '.', 0)

    # =============================================
    # Matrix and matrix metadate output to files
    # =============================================
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

    return matrix, queries, targets, matrix_file, metadata_file

def load_alignment_matrix_from_file(matrix_path: str,
                                    metadata_path: str):
    """
    Loads alignment matrix and matrix metadata into memory

    Args:
        matrix_path (str): Path to .npy file containing alignment matrix
        metadata_path (str): Path to matrix metadata file

    Returns:
        matrix (numpy array): Alignment matrix
        queries (list): List of query protein IDs
        targets (list): List of reference protein IDs
    """
    matrix = np.load(matrix_path)

    with open(metadata_path, "r") as file:
        metadata = json.load(file)

    queries = metadata["queries"]
    targets = metadata["targets"]

    logger.info(f"Loaded matrix shape: {matrix.shape}")
    logger.info(f"Loaded {len(queries)} queries and {len(targets)} targets")

    return matrix, queries, targets

def create_embedding_file(output_file: str,
                          embedding: np.ndarray,
                          queries: list,
                          clusters: np.ndarray,
                          col_names: list,
                          ):
    """
    Creates a TSV file containing protein embeddings with metadata and cluster assignments.

    Args:
        output_file (str): Path to the output TSV file to be created
        embedding (np.ndarray): 2D array of embedding vectors, shape (n_proteins, n_dimensions)
        queries (list): List of protein IDs in format "genomeID_locusNr" (e.g., "genome_ABC_001")
        clusters (np.ndarray): 1D array of cluster assignments for each protein
        col_names (list): Column names for each embedding dimension (e.g., ["dim_0", "dim_1", ...])
        metadata_protein (str): Protein metadata
        metadata_genome (str): Genome metadata
    """
    logger.info(f"Creating embedding TSV with {len(queries)} proteins")

    loci = [seqID.rsplit("_", 1)[1] for seqID in queries]
    genome_ids = [seqID.rsplit("_", 1)[0] for seqID in queries]

    with open(output_file, "w", encoding="utf-8") as f:
        base_cols = ["seqID", "locus_nr", "genome_ID"] + col_names + ["cluster"]

        f.write("\t".join(base_cols) + "\n")
        for q, c, l, g, emb in zip(queries, clusters, loci, genome_ids, embedding):
            row = [q, l, g] + [f"{x:.6f}" for x in emb] + [str(c)]
            f.write("\t".join(row) + "\n")

    logger.info(f"Finished writing TSV: {output_file}")


def create_embedding_dataframe(embedding: np.ndarray,
                               queries: list,
                               clusters: np.ndarray,
                               col_names: list
                               ):
    """
    Creates DataFrame from t-SNE embedding results with clustering and optional metadata.

    Args:
        embedding (np.ndarray): t-SNE embedding coordinates
        queries (list): List of protein IDs corresponding to embedding rows
        clusters (np.ndarray): Cluster assignments for each protein
        col_names (list): Column names for embedding dimensions
        metadata_protein (str, optional): Path to protein metadata TSV file
        metadata_genome (str, optional): Path to genome metadata TSV file

    Returns:
        df (pd.DataFrame): DataFrame with protein IDs, coordinates, clusters, and metadata
    """
    logger.debug(f"Creating embedding dataframe with {len(queries)} proteins")
    df = pd.DataFrame(embedding, columns=col_names)
    df['seqID'] = queries
    df['cluster'] = clusters

    logger.debug("Extracting locus and genome information from protein IDs")

    # extract locus and genome info from protein IDs
    df['locus_nr'] = df['seqID'].astype(str).str.rsplit("_", n=1).str[1]
    df['genome_ID'] = df['seqID'].astype(str).str.rsplit("_", n=1).str[0]

    # reorder columns
    base_cols = ['seqID', 'locus_nr', 'genome_ID'] + col_names + ['cluster']
    df = df[base_cols]

    logger.debug(f"Final dataframe shape: {df.shape}")
    return df


def tsne_embedding(matrix: np.ndarray,
                   queries: list,
                   basename: str,
                   perplexity: int,
                   iterations: int,
                   exaggeration: int,
                   threads: int,
                   force: bool = False
                   ):
    """
    Performs t-SNE embedding and DBSCAN clustering on alignment matrix.

    Creates early (exaggerated) and final t-SNE embeddings with DBSCAN clustering
    applied to the early embedding. Same cluster labels are used for both embeddings.

    Args:
        matrix (np.ndarray): Alignment matrix for embedding
        queries (list): List of protein IDs corresponding to matrix rows
        basename (str): Base name for output files
        perplexity (int): t-SNE perplexity parameter
        iterations (int): Number of optimization iterations per phase
        exaggeration (int): Early exaggeration factor
        threads (int): Number of threads for parallel processing
        force (bool): Overwrite existing files if True

    Returns:
        early_filename: Path to early embedding TSV file
        final_filename: Path to final embedding TSV file
    """
    logger.info("Starting t-SNE embedding process")
    logger.info(f"Matrix shape: {matrix.shape}")
    logger.info(f"Perplexity: {perplexity}, Iterations: {iterations}, Exaggeration: {exaggeration}")
    logger.info(f"Using {threads} threads")

    # ========================
    # Multiscale affinites
    # ========================

    # compute affinity matrix using multiple perplexity values for local and global structure
    logger.info("Computing multiscale affinities")
    affinities = openTSNE.affinity.Multiscale(matrix,
                                              perplexities=[20, perplexity],
                                              metric="cosine",
                                              n_jobs=threads)

    # ========================
    # Initialize embedding
    # ========================

    # initialize with PCA for better convergence
    logger.info("Initializing t-SNE with PCA")
    pca_init = openTSNE.initialization.pca(matrix)

    # create t-SNE embedding with auto learning rate and single DOF
    tsne_embed = openTSNE.TSNEEmbedding(pca_init,
                                        affinities,
                                        n_jobs=threads,
                                        dof=1,
                                        learning_rate="auto")


    # ============================
    # Early optimization phase
    # ============================

    # optimize with exaggeration to spread clusters; lower momentum for exploration
    logger.info(f"Starting early optimization with exaggeration={exaggeration}")
    tsne_embed.optimize(n_iter=iterations,
                        exaggeration=exaggeration,
                        momentum=0.5,
                        inplace=True,
                        n_jobs=threads)
    logger.info("Early optimization completed")

    # =================================
    # Clustering on early embedding
    # =================================

    # apply DBSCAN to identify protein clusters
    logger.info("Performing DBSCAN clustering on early embedding")
    clustering = DBSCAN(eps=1, min_samples=10).fit(tsne_embed)

    # calculate cluster statistics (label -1 indicates noise points)
    n_clusters = len(set(clustering.labels_)) - (1 if -1 in clustering.labels_ else 0)
    n_noise = list(clustering.labels_).count(-1)
    logger.info(f"DBSCAN clustering results: {n_clusters} clusters, {n_noise} noise points")

    # ================================
    # Save early embedding results
    # ================================
    logger.debug("Saving early embedding results")
    early_filename = ensure_path(target=f"{basename}_tsne_early_clust.tsv", force=force)
    create_embedding_file(early_filename, tsne_embed, queries, clustering.labels_, ['tsne1', 'tsne2'])
    logger.info(f"Early embedding saved to: {early_filename}")

    # ========================
    # Final optimization phase
    # ========================

    # fine-tune embedding without exaggeration; higher momentum for stability
    logger.info("Starting final optimization")
    tsne_embed.optimize(n_iter=iterations,
                        momentum=0.8,
                        inplace=True,
                        n_jobs=threads)
    logger.info("Final optimization completed")

    # ================================
    # Save final embedding results
    # ===============================
    logger.debug("Saving final embedding results")
    final_filename = ensure_path(target=f"{basename}_tsne_final_clust.tsv", force=force)
    create_embedding_file(final_filename, tsne_embed, queries, clustering.labels_, ['tsne1', 'tsne2'])
    logger.info(f"Final embedding saved to: {final_filename}")

    return early_filename, final_filename

def fetch_protein_metadata(db_path: str,
                           protein_ids: list,
                           column: str,
                           batch_size: int = 500
                           ):
    valid_cols = ['COG_ID', 'KEGG_ID', 'Pfam_ID']

    if column not in valid_cols:
        logger.warning(f"No valid protein metadata columns requested")
        return pd.DataFrame({'seqID': protein_ids})

    conn = sqlite3.connect(db_path)

    col_str = ', '.join(['seqID', column])

    all_results = []
    n_batches = (len(protein_ids) + batch_size - 1) // batch_size
    logger.info(f"Fetching protein metadata in {n_batches} batches")

    for i in range(0, len(protein_ids), batch_size):
        batch = protein_ids[i:i + batch_size]
        placeholders = ','.join(['?'] * len(batch))
        query = f"""SELECT {col_str} from all_data where seqID in ({placeholders})"""

        batch_df = pd.read_sql_query(query, conn, params=batch)
        all_results.append(batch_df)

    conn.close()

    df = pd.concat(all_results, ignore_index=True)

    logger.info(f"Retrieved {len(df)} genome metadata records total")

    return df


def fetch_genome_metadata(db_path: str,
                          genome_ids: list,
                          column: str,
                          batch_size: int = 500
                           ):
    valid_cols = ['domain', 'phylum', 'class', 'order_tax', 'family', 'genus', 'species']

    if column not in valid_cols:
        logger.warning(f"No valid protein metadata columns requested")
        return pd.DataFrame({'genome_ID': genome_ids})

    conn = sqlite3.connect(db_path)
    col_str = ', '.join(['genome_ID',  column])

    # Process in batches
    all_results = []
    n_batches = (len(genome_ids) + batch_size - 1) // batch_size
    logger.info(f"Fetching genome metadata in {n_batches} batches")

    for i in range(0, len(genome_ids), batch_size):
        batch = genome_ids[i:i + batch_size]
        placeholders = ','.join(['?'] * len(batch))
        query = f"SELECT {col_str} FROM taxonomy WHERE genome_ID IN ({placeholders})"

        batch_df = pd.read_sql_query(query, conn, params=batch)
        all_results.append(batch_df)
        logger.debug(f"Batch {len(all_results)}/{n_batches}: retrieved {len(batch_df)} records")

    conn.close()

    df = pd.concat(all_results, ignore_index=True) if all_results else pd.DataFrame()

    logger.info(f"Retrieved {len(df)} genome metadata records total")
    return df


def extend_embedding_with_metadata(df: pd.DataFrame,
                                   db_path: str,
                                   protein_col: str = None,
                                   genome_col: str = None):
    if not db_path:
        logger.warning("No database path provided, skipping metadata enrichment")
        return df

    if protein_col is not None:
        logger.info(f"Fetching protein metadata: {protein_col}")
        protein_metadata = fetch_protein_metadata(db_path, df['seqID'].tolist(), protein_col)
        df = pd.merge(df, protein_metadata, on='seqID', how='left')

    if genome_col is not None:
        logger.info(f"Fetching genome metadata: {genome_col}")
        genome_metadata = fetch_genome_metadata(db_path, df['genome_ID'].unique().tolist(), genome_col)
        df = pd.merge(df, genome_metadata, on='genome_ID', how='left')

    return df


def plot_clusters(tsv_file: str,
                  output: str,
                  db_path: str = None,
                  metadata_protein: str = None,
                  metadata_genome: str = None,
                  svg: bool = False,
                  force: bool = False,
                  show_cluster_numbers: bool = False
                  ):
    """
    Create a t-SNE scatter plot visualization of DBSCAN clustering results.

    Reads clustering results from a TSV file and generates a scatter plot showing
    the t-SNE embedding with different colors for each cluster. Noise points
    (cluster -1) are displayed in light gray with reduced opacity.

    Args:
        tsv_file (str): Path to input TSV file containing tSNE coordinates and cluster assignment
        output (str): Base path for output plot file (adds "_tsne_clusters" suffix)
        svg (bool): Generate plot in SVG format
        force (bool): Overwrite existing files if True
        show_cluster_numbers (bool): Display cluster number on cluster centers in output plot
    """
    # load clustering results from TSV file
    logger.info(f"Creating t-SNE plot from: {tsv_file}")
    df = pd.read_csv(tsv_file, sep='\t')
    logger.info(f"Loaded {len(df)} data points for plotting")

    if 'seqID' not in df.columns:
        raise KeyError("Expected column 'seqID' not found in TSV")

    df['seqID'] = df['seqID'].astype(str).str.strip()
    df['genome_ID'] = (
        df['seqID']
        .apply(parse_protein_identifier)
        .astype(str)
        .str.strip()
    )

    color_column = None
    if db_path:
        if metadata_protein:
            logger.info(f"Fetching protein metadata column '{metadata_protein}' for coloring")
            df = extend_embedding_with_metadata(df, db_path, protein_col=metadata_protein, genome_col=None)
            color_column = metadata_protein
        elif metadata_genome:
            logger.info(f"Fetching genome metadata column '{metadata_genome}' for coloring")
            df = extend_embedding_with_metadata(df, db_path, protein_col=None, genome_col=metadata_genome)
            color_column = metadata_genome

    if (metadata_protein or metadata_genome) and not db_path:
        logger.info("Metadata selected but path to SQLite DB not provided, coloring by cluster")
        color_column = 'cluster'
    elif color_column is None or color_column not in df.columns:
        logger.info("No metadata column selected, coloring by cluster")
        color_column = 'cluster'

    # ========================
    # Initialize plot
    # ========================
    plt.figure(figsize=(14, 10))

    if color_column == 'cluster':
        # create masks to separate noise points (cluster -1) from clustered points
        noise_mask = df['cluster'] == -1
        clustered_mask = df['cluster'] != -1

        n_noise = noise_mask.sum()
        n_clustered = clustered_mask.sum()
        logger.info(f"Plotting {n_clustered} clustered points and {n_noise} noise points")

        # ========================
        # Plot noise points
        # ========================

        # display unclustered points in light gray with low opacity
        if noise_mask.any():
            plt.scatter(df.loc[noise_mask, 'tsne1'],
                        df.loc[noise_mask, 'tsne2'],
                        c='lightgray',
                        alpha=0.4,
                        s=12,
                        edgecolors='none',
                        label='Noise')

        # ========================
        # Plot clustered points
        # ========================
        if clustered_mask.any():
            unique_clusters = df.loc[clustered_mask, 'cluster'].unique()

            # use tab20 colormap for up to 20 clusters, switch to hsv for more
            colors = plt.cm.tab20(np.linspace(0, 1, len(unique_clusters)))
            if len(unique_clusters) > 20:
                colors = plt.cm.hsv(np.linspace(0, 1, len(unique_clusters)))

            logger.info(f"Plotting {len(unique_clusters)} unique clusters")

            # plot each cluster with unique color
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

                # optionally annotate cluster centers with cluster numbers
                if show_cluster_numbers:
                    centroid_x = df.loc[cluster_mask, 'tsne1'].mean()
                    centroid_y = df.loc[cluster_mask, 'tsne2'].mean()
                    plt.annotate(str(cluster),
                                xy=(centroid_x, centroid_y),
                                fontsize=16,
                                fontweight='bold',
                                ha='center',
                                va='center',
                                color='black',
                                alpha=1)
    else:
        unique_vals = df[color_column].dropna().unique()
        n_vals = len(unique_vals)

        na_mask = df[color_column].isna()
        if na_mask.any():
            plt.scatter(df.loc[na_mask, 'tsne1'], df.loc[na_mask, 'tsne2'],
                        c='lightgray', alpha=0.4, s=12, edgecolors='none', label='No data')

        if df[color_column].dtype in ['object', 'category']:
            colors = plt.cm.tab20(np.linspace(0, 1, n_vals))
            if n_vals > 20:
                colors = plt.cm.hsv(np.linspace(0, 1, n_vals))
            for i, val in enumerate(unique_vals):
                mask = df[color_column] == val
                plt.scatter(df.loc[mask, 'tsne1'], df.loc[mask, 'tsne2'],
                            c=[colors[i]], s=15, alpha=0.7,
                            edgecolors='white', linewidths=0.3, label=str(val))
        else:
            scatter = plt.scatter(df['tsne1'], df['tsne2'],
                                  c=df[color_column], cmap='viridis', s=15,
                                  alpha=0.7, edgecolors='white', linewidths=0.3)
            plt.colorbar(scatter, label=color_column)

    # =============================
    # Configure plot appearance
    # =============================
    plt.xlabel('t-SNE 1', fontsize=12)
    plt.ylabel('t-SNE 2', fontsize=12)
    plt.title('t-SNE Embedding with DBSCAN Clusters', fontsize=14)

    if color_column == 'cluster' and 'cluster' in df.columns:
        unique_clusters = df[df['cluster'] != -1]['cluster'].unique()
        if len(unique_clusters) <= 60:
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
        else:
            logger.info("Cluster number exceeds 60, omitting legend.")
    elif color_column != 'cluster':
        unique_vals = df[color_column].dropna().unique()
        if len(unique_vals) <= 60:
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
        else:
            logger.info("Cluster number exceeds 60, omitting legend.")

    plt.tight_layout()


    # =============================
    # Determine output filename
    # =============================

    # parse input filename to determine if early or final embedding
    if 'early_clust' in tsv_file:
        prefix = determine_dataset_name(tsv_file, '.', 0, '_tsne_early_clust')
        base_name = f"{prefix}_tsne_early"
    elif 'final_clust' in tsv_file:
        prefix = determine_dataset_name(tsv_file, '.', 0, '_tsne_final_clust')
        base_name = f"{prefix}_tsne_final"
    else:
        raise ValueError(f"Unexpected TSV filename: {tsv_file}")

    if color_column and color_column != 'cluster':
        filename = f"{base_name}_by_{color_column}.{ 'svg' if svg else 'png' }"
    else:
        filename = f"{base_name}_clusters.{ 'svg' if svg else 'png' }"


    # ========================
    # Save plot
    # ========================
    output_file = ensure_path(output, filename, force=force)
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    logger.info(f"Plot saved to: {output_file}")


def pick(final_embedding_file: str,
         fasta: str,
         no_cluster: int,
         output: str,
         force: bool = False):
    """
    Extract proteins belonging to a specific cluster into a separate FASTA file.

    Reads the t-SNE clustering results and extracts all protein sequences that
    belong to the specified cluster number, writing them to a new FASTA file.

    Args:
        final_embedding_file (str): Path to TSV file with t-SNE embeddings and cluster assignments
        fasta (str): Path to input FASTA file containing all protein sequences
        no_cluster (int): Cluster number to extract
        output (str): Output directory path
        force (bool): Overwrite existing files if True

    Returns:
        str: Path to the output FASTA file containing selected cluster sequences
    """
    # ===========================
    # Prepare output filename
    # ===========================

    # generate output filename based on input file and cluster number
    prefix = determine_dataset_name(final_embedding_file, '.', 0, '_tsne_final_clust')
    cluster_fasta = ensure_path(output, f"{prefix}_cluster_{no_cluster}.faa", force=force)

    # ============================================
    # Extract protein IDs from target cluster
    # ============================================
    cluster_seqIDs = []
    with open(final_embedding_file, 'r') as file:
        # parse header to find column indices
        header = file.readline().strip().split()
        cluster_idx = header.index("cluster")
        seqIDx = header.index("seqID")

        # iterate through data rows and collect protein IDs matching target cluster
        for line in file:
            parts = line.strip().split('\t')
            seqID = parts[seqIDx]
            prot_cluster_no = parts[cluster_idx]

            if int(prot_cluster_no) == int(no_cluster):
                cluster_seqIDs.append(seqID)

    # ==================================
    # Load original FASTA sequences
    # ==================================
    input_fasta_dict = read_fasta_to_dict(fasta)

    # ============================================
    # Write cluster sequences to output FASTA
    # ============================================

    # write sequences for all proteins in the selected cluster
    with open(cluster_fasta, 'w') as file:
        for id in cluster_seqIDs:
            if id in input_fasta_dict:
                file.write(f">{id}\n{input_fasta_dict[id]}\n")

    return cluster_fasta

# ==================================================
#
#                   CLI FUNCTIONS
#
# ==================================================

def matrix(fasta: str,
          output: str,
          subset: str,
          subset_size: int,
          threads: int,
          force: bool = False
          ):
    """
    Build alignment matrix from FASTA sequences using DIAMOND alignment.

    Creates a subset of sequences (either provided or randomly sampled), performs
    DIAMOND alignment against the full dataset, and constructs an alignment matrix
    from the results. The matrix represents alignment scores between sequences.

    Args:
        fasta (str): Path to input FASTA file with protein sequences
        output (str): Base name for output files
        subset (str): Path to subset FASTA file (if None, creates random subset)
        subset_size (int): Number of sequences for subset (ignored if subset provided)
        threads (int): Number of threads for DIAMOND alignment
        force (bool): Overwrite existing files if True

    Returns:
        matrix_file (str): Path to .npy file containing numpy matrix
        metadata_file (str): Path to .json file containing matrix metadata (queries, targets, matrix stats)
    """
    logger.info("=== Starting Matrix Construction ===")
    logger.info(f"Input FASTA: {fasta}")
    logger.info(f"Output basename: {output}")
    logger.info(f"Threads: {threads}")
    if subset:
       logger.info(f"Using provided subset: {subset}")
       subset_fasta = subset
       subset_size = count_fasta_sequences(subset_fasta)
    else:
       logger.info("Creating random subset from input FASTA")
       subset_fasta = fasta_subsample(fasta, output, subset_size, force=force)

    logger.info("=== Phase 1: DIAMOND Alignment ===")
    align_output = run_diamond_alignment(fasta, subset_fasta, subset_size, threads, force=force)

    logger.info("=== Phase 2: Matrix Construction ===")
    _, _, _, matrix_file, metadata_file = build_alignment_matrix_split(align_output, output, force)

    return matrix_file, metadata_file


def cluster(matrix_path: str,
         matrix_metadata_path: str,
         output: str,
         perplexity: int = 50,
         iterations: int = 500,
         exaggeration: int = 6,
         threads: int = 1,
         force: bool = False,
         ):
    """
    Perform t-SNE embedding and DBSCAN clustering on alignment matrix.

    Creates early (exaggerated) and final t-SNE embeddings with DBSCAN clustering
    applied to the early embedding. Same cluster labels are used for both embeddings.
    Optionally enriches results with protein and genome metadata.

    Args:
        matrix_path (str): Path to alignment matrix file
        matrix_metadata_path (str): Path to matrix metadata file
        output (str): Base name for output files
        perplexity (int): t-SNE perplexity parameter
        iterations (int): Number of optimization iterations per phase
        exaggeration (int): Early exaggeration factor
        threads (int): Number of threads for parallel processing
        force (bool): Overwrite existing files if True

    Returns:
        early_filename (str), final_filename (str): Paths to early and final embedding TSV files
    """
    logger.info("=== Starting t-SNE Embedding ===")
    logger.info(f"Matrix file: {matrix_path}")
    logger.info(f"Metadata file: {matrix_metadata_path}")
    if output:
       logger.info(f"Output basename: {output}")

    prefix = determine_dataset_name(matrix_path, '.', 0, '_matrix')

    matrix, queries, targets = load_alignment_matrix_from_file(matrix_path, matrix_metadata_path)

    early_filename, final_filename = tsne_embedding(matrix=matrix,
                              queries=queries,
                              basename=prefix,
                              perplexity=perplexity,
                              iterations=iterations,
                              exaggeration=exaggeration,
                              threads=threads,
                              force=force,
                              )
    logger.info("=== t-SNE Embedding Completed ===")
    return early_filename, final_filename


def casm_plot(early_clust_path: str,
        full_clust_path: str,
        output: str,
        db_path: str = None,
        metadata_protein: str = None,
        metadata_genome: str = None,
        svg: bool = False,
        force: bool = False,
        show_cluster_numbers: bool = False):
    """
    Generate t-SNE cluster visualization plots for early and final embeddings.

    Creates scatter plots showing DBSCAN clustering results for both early
    (exaggerated) and final t-SNE embeddings. Generates two plot files with
    "_early_tsne_clusters" and "_final_tsne_clusters" suffixes.

    Args:
        early_clust_path (str): Path to early clustering results TSV file
        full_clust_path (str): Path to final clustering results TSV file
        output (str): Base name for output plot files

    Returns:
        None: Saves plots to disk as files
    """
    logger.info("=== Starting Plot Generation ===")
    logger.info(f"Early clustering file: {early_clust_path}")
    logger.info(f"Full clustering file: {full_clust_path}")
    logger.info(f"Output basename: {output}")

    # Generate plots
    logger.info("Generating early clustering plot")
    plot_clusters(early_clust_path,
                  output=output,
                  db_path=db_path,
                  metadata_protein=metadata_protein,
                  metadata_genome=metadata_genome,
                  force=force,
                  svg=svg,
                  show_cluster_numbers=show_cluster_numbers)

    logger.info("Generating full clustering plot")
    plot_clusters(full_clust_path,
                  output=output,
                  db_path=db_path,
                  metadata_protein=metadata_protein,
                  metadata_genome=metadata_genome,
                  force=force,
                  svg=svg,
                  show_cluster_numbers=show_cluster_numbers)

    logger.info("=== Plot Generation Completed ===")


def casm(fasta: str,
         output: str,
         subset: str,
         subset_size: int,
         db_path: str = None,
         threads: int = 1,
         perplexity: int = 50,
         iterations: int = 500,
         exaggeration: int = 6,
         metadata_protein: str = None,
         metadata_genome: str = None,
         svg: bool = False,
         force: bool = False,
         show_cluster_numbers: bool = False
         ):
    """
    Run complete CASM analysis pipeline.

    Args:
        fasta (str): Path to input FASTA file
        output (str): Output basename
        subset (str, optional): Path to subset FASTA file. If None, creates random subset
        subset_size (int): Size of subset for alignment
        threads (int): Number of threads
        perplexity (int): t-SNE perplexity parameter
        iterations (int): Number of optimization iterations
        exaggeration (int): Early exaggeration parameter
        metadata_protein (str, optional): Path to protein metadata file
        metadata_genome (str, optional): Path to genome metadata file
        svg (bool): Generate plot in SVG format
        force (bool): Force overwrite existing files
        show_cluster_numbers (bool): Display cluster enumeration in cluster centers in output plot

    Returns:
        sum_dict: Dictionary containing paths to all generated files
    """
    logger.info("=== Starting Complete CASM Analysis ===")
    logger.info(f"Input FASTA: {fasta}")
    logger.info(f"Output basename: {output}")
    logger.info(f"Subset size: {subset_size}")
    logger.info(f"Threads: {threads}")

    # Phase 1: Matrix construction
    logger.info("=== Phase 1: Matrix Construction ===")
    matrix_file, metadata_file = matrix(fasta=fasta,
                                        output=output,
                                        subset=subset,
                                        subset_size=subset_size,
                                        threads=threads,
                                        force=force
                                        )

    # Phase 2: t-SNE embedding
    logger.info("=== Phase 2: t-SNE Embedding ===")
    early_filename, final_filename = cluster(
        matrix_path=matrix_file,
        matrix_metadata_path=metadata_file,
        output=output,
        perplexity=perplexity,
        iterations=iterations,
        exaggeration=exaggeration,
        threads=threads,
        force=force
    )

    # Phase 3: Plot generation
    logger.info("=== Phase 3: Plot Generation ===")
    casm_plot(
        early_clust_path=early_filename,
        full_clust_path=final_filename,
        output=output,
        db_path=db_path,
        metadata_protein=metadata_protein,
        metadata_genome=metadata_genome,
        svg=svg,
        force=force,
        show_cluster_numbers=show_cluster_numbers
    )

    logger.info("=== Complete CASM Analysis Completed ===")

    sum_dict = {
        "matrix_file": matrix_file,
        "metadata_file": metadata_file,
        "early_embedding": early_filename,
        "final_embedding": final_filename
    }

    return sum_dict



