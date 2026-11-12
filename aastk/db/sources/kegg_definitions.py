import logging
logger = logging.getLogger(__name__)

def parse_kegg_definitions_file(kegg_def_filepath: str) -> list:
    kegg_definitions_data = []
    with open(kegg_def_filepath, 'r') as f:
        f.readline()
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 12:
                continue

            kegg_id = parts[0]
            kegg_definition = parts[11]

            kegg_definitions_data.append((kegg_id, kegg_definition))

    return kegg_definitions_data
