# input: protein fasta
# optional arguments: output dir, model, batch size (sequences) depending on sequence length?, force/keep
# output: .h5 (real embedding output), metadata json?, tsv for quick inspection?

from pathlib import Path

#import torch
#import h5py
#from transformers import AutoTokenizer, T5EncoderModel

from aastk.util import ensure_path

import logging
logger = logging.getLogger(__name__)

MODEL_REGISTRY = {
    "prott5": "Rostlab/prot_t5_xl_half_uniref50-enc"
}

MAX_RESIDUES_PER_BATCH = 5000

def resolve_output_path(
        fasta_path: Path,
        output: str = None,
        force: bool = False):

    fasta_name = fasta_path.name

    if fasta_name.endswith(".gz"):
        fasta_name = Path(fasta_name).stem

    fasta_stem = Path(fasta_name).stem
    output_file = f"{fasta_stem}_embeddings.h5"

    # no output provided as argument
    if output is None:
        return Path(ensure_path(
            path=str(fasta_path.parent),
            target=output_file,
            force=force
        ))

    output_path = Path(output)

    # output file name is provided
    if output_path.suffix.lower() in [".h5", ".hdf5"]:
        return Path(ensure_path(
            path=str(output_path.parent),
            target=output_path.name,
            force=force
        ))

    # output is something else -> handled as directory
    return Path(ensure_path(
            path=str(output_path),
            target=output_file,
            force=force
    ))

def collect_fasta_records(fasta_path: Path):
    """
    reading protein FASTA records and preparing them
    """

def load_plm_model():
    """
    loading PLM model and respective Tokenizer
    """

def write_embeddings():
    """
    writing mean-pooled protein embeddings to h5
    repeatedly call embed_current_batch
    """

def batch_handling():
    """
    creating batches while limiting sequence/residue count
    """

def embed_current_batch():
    """
    do model call for current batch and return embeddings
    """


def plm_embedder(
        fasta: str,
        model: str,
        output: str = None,
        batch_size: int = 1,
        force: bool = False):
    logger.info("PLM embedder reached successfully")

    fasta_path = Path(fasta)
    output_path = resolve_output_path(fasta_path=fasta_path, output = output)

