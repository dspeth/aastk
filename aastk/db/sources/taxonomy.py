import gzip
from pathlib import Path

from aastk.db.schema import TAXONOMY_COLUMNS

import logging
logger = logging.getLogger(__name__)

def parse_taxonomy_file(taxonomy_filepath: str) -> list:
    filepath = Path(taxonomy_filepath)

    if filepath.suffix == '.gz':
        file_handle = gzip.open(filepath, 'rt')
    else:
        file_handle = open(filepath, 'r')

    taxonomy_data = []
    with file_handle as f:
        next(f)

        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < len(TAXONOMY_COLUMNS):
                continue

            taxonomy_data.append(tuple(parts[i] for i in range(len(TAXONOMY_COLUMNS))))

    return taxonomy_data
