# prepares input from an hdf5/h5 file to work with the CASM workflow
# imports all reusable downstream functions from the original casm.py

# ATTENTION / notes for later
# protein IDs used in CASM might not be identical to the input IDs (if they contained "/")
# load function assumes every root is an embedding (pay attention if we add e.g. metadata)
# all embeddings are loaded into memory

import logging
from pathlib import Path

import h5py
import numpy as np

from .casm import *

logger = logging.getLogger(__name__)

def load_plm_embeddings_from_h5(embedding_path: str):
    """
    load mean-pooled embeddings from an hdf5 file

    returns:
        matrix (numpy array)
        protein IDs (list): corresponding to the matrix rows
    """
    with h5py.File(embedding_path, "r") as h5_file:
        protein_ids = sorted(h5_file.keys())

        if not protein_ids:
            raise ValueError(
                f"No embeddings found in HDF5 file: {embedding_path}"
            )

        embeddings = []

        for protein_id in protein_ids:
            embedding = np.asarray(
                h5_file[protein_id][()],
                dtype=np.float32
            )
            embeddings.append(embedding)

    matrix = np.stack(embeddings, axis=0)

    logger.info(f"Loaded {len(protein_ids)} protein embeddings")
    logger.info(f"Embedding matrix shape: {matrix.shape}")

    return matrix, protein_ids

def casm_plm(
        embeddings: str,
        output: str,
        db_path: str = None,
        threads: int = 1,
        perplexity: int = 50,
        iterations: int = 500,
        exaggeration: int = 6,
        metadata: str = None,
        keep: bool = False,
        svg: bool = False,
        force: bool = False,
        show_cluster_numbers: bool = False):
    """
    returns dictionary containing all generated file paths
    """
    logger.info("=== starting PLM CASM analysis ===")
    logger.info(f"embedding file: {embeddings}")

    matrix, protein_ids = load_plm_embeddings_from_h5(embedding_path=embeddings)

    basename = Path(embeddings).stem.removesuffix("_embeddings")

    early_filename, final_filename = tsne_embedding(
        matrix=matrix,
        queries=protein_ids,
        output=output,
        basename=basename,
        perplexity=perplexity,
        iterations=iterations,
        exaggeration=exaggeration,
        threads=threads,
        force=force
    )

    early_cluster_plot = casm_plot(
        clust_path=early_filename,
        output=output,
        metadata=metadata,
        db_path=db_path,
        svg=svg,
        force=force,
        show_cluster_numbers=show_cluster_numbers
    )

    full_cluster_plot = casm_plot(
        clust_path=final_filename,
        output=output,
        metadata=metadata,
        db_path=db_path,
        svg=svg,
        force=force,
        show_cluster_numbers=show_cluster_numbers
    )

    logger.info("=== PLM CASM analysis completed ===")

    return {
        "early_filename": early_filename,
        "final_filename": final_filename,
        "early_cluster_plot": early_cluster_plot,
        "final_cluster_plot": full_cluster_plot
    }