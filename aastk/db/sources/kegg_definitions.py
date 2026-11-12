import logging

from aastk.db.sources import QUOTE_TRANSLATION, read_lines_utf8_or_cp1252

logger = logging.getLogger(__name__)

def parse_kegg_definitions_file(kegg_def_filepath: str) -> list:
    kegg_definitions_data = []
    for line in read_lines_utf8_or_cp1252(kegg_def_filepath)[1:]:
        parts = line.split('\t')
        if len(parts) < 12:
            continue

        kegg_id = parts[0]
        kegg_definition = parts[11].translate(QUOTE_TRANSLATION)

        kegg_definitions_data.append((kegg_id, kegg_definition))

    return kegg_definitions_data
