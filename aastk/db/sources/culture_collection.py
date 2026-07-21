import logging
logger = logging.getLogger(__name__)

def parse_culture_collection_file(cc_filepath: str) -> list:
    culture_collection_data = []
    with open(cc_filepath, 'r') as f:
        next(f)

        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue

            genome_id = parts[0]
            culture_collection = int(parts[1])

            culture_collection_data.append((genome_id, culture_collection))

    return culture_collection_data
