import logging

logger = logging.getLogger(__name__)

def annotate(
        db_fasta: str,
        yaml: str,
        output: str,
        query_protein: str | None = None,
        query_genome: str | None = None,
        keep: bool = False,
        force: bool = False):
    logger.info("connection to annotate worked")