import logging
logger = logging.getLogger(__name__)

def parse_cog_definitions_file(cog_def_filepath: str) -> list:
    cog_definitions_data = []
    with open(cog_def_filepath, 'r') as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue

            cog_id = parts[0]
            cog_definition = parts[2]

            cog_definitions_data.append((cog_id, cog_definition))

    return cog_definitions_data
