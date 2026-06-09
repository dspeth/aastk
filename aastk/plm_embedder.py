# input: protein fasta
# optional arguments: output dir, model, batch size (sequences) depending on sequence length?, force/keep
# output: .h5 (real embedding output), metadata json?, tsv for quick inspection?

import logging
logger = logging.getLogger(__name__)

def plm_embedder(
        fasta: str,
        model: str,
        output: str = None,
        batch_size: int = 1,
        force: bool = False):
    logger.info("PLM embedder reached successfully")
