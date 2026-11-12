import logging

from aastk.db.sources import QUOTE_TRANSLATION, read_lines_utf8_or_cp1252

logger = logging.getLogger(__name__)

def parse_cog_definitions_file(cog_def_filepath: str) -> list:
    cog_definitions_data = []
    for line in read_lines_utf8_or_cp1252(cog_def_filepath):
        parts = line.split('\t')
        if len(parts) < 3:
            continue

        cog_id = parts[0]
        cog_definition = parts[2].translate(QUOTE_TRANSLATION)

        cog_definitions_data.append((cog_id, cog_definition))

    return cog_definitions_data
