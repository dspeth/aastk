import csv
import logging
from pathlib import Path

import yaml

from aastk.util import ensure_path

logger = logging.getLogger(__name__)

POSITION_INDICES = list(range(-7, 0)) + list(range(1, 10))
PUBLICATION_INDICES = range(1, 6)


def _get(row: dict, column: str) -> str:
    value = row.get(column, 'NA')
    return value.strip() if value and value.strip() else 'NA'


def _as_float(value: str):
    return float(value) if value != 'NA' else 'NA'


def _as_int(value: str):
    return int(value) if value != 'NA' else 'NA'


def build_info_yaml(row: dict, output_dir: str, force: bool = False) -> str:
    gdp_id = row['GDP ID']

    data = {
        gdp_id: {
            'gene_family': _get(row, 'gene family'),
            'description': _get(row, 'description'),
            'alternative_names': _get(row, 'alternative_names'),
            'database_citation': _get(row, 'database_citation'),
            'version': {
                'globdb': _as_int(_get(row, 'globdb_version')),
                'gene_family': _as_int(_get(row, 'gene_family_version')),
            },
            'annotation': {
                'COG': _get(row, 'COG_annotation'),
                'KEGG': _get(row, 'KEGG_annotation'),
                'PFAM': _get(row, 'PFAM_annotation'),
            },
            'cutoffs': {
                'lasr': _as_float(_get(row, 'lasr_cutoff')),
                'selfmax': _as_int(_get(row, 'selfmax_cutoff')),
                'selfmin': _as_int(_get(row, 'selfmin_cutoff')),
                'matrix': _get(row, 'matrix'),
            },
        }
    }

    yaml_path = ensure_path(output_dir, f'{gdp_id}/info.yaml', force=force)
    with open(yaml_path, 'w') as f:
        yaml.dump(data, f, sort_keys=False, default_flow_style=False)

    return yaml_path


def build_synteny_yaml(row: dict, output_dir: str, force: bool = False):
    gdp_id = row['GDP ID']
    positions = {}

    for position in POSITION_INDICES:
        gene_id = _get(row, f'pos_{position}_ID')
        if gene_id == 'NA':
            continue
        fraction = _as_float(_get(row, f'pos_{position}_fraction'))
        positions[position] = [{gene_id: fraction}]

    if not positions:
        return None

    data = {gdp_id: positions}

    yaml_path = ensure_path(output_dir, f'{gdp_id}/synteny.yaml', force=force)
    with open(yaml_path, 'w') as f:
        yaml.dump(data, f, sort_keys=False, default_flow_style=False)

    return yaml_path


def build_validation_yaml(row: dict, output_dir: str, force: bool = False):
    gdp_id = row['GDP ID']
    publications = {}

    for i in PUBLICATION_INDICES:
        pub_id = _get(row, f'pub_{i}_ID')
        if pub_id == 'NA':
            continue

        sequences = _get(row, f'pub_{i}_sequences')
        publications[pub_id] = {
            'title': _get(row, f'pub_{i}_title'),
            'doi': _get(row, f'pub_{i}_doi'),
            'sequences': [s.strip() for s in sequences.split(',')] if sequences != 'NA' else 'NA',
        }

    if not publications:
        return None

    data = {gdp_id: publications}

    yaml_path = ensure_path(output_dir, f'{gdp_id}/validation.yaml', force=force)
    with open(yaml_path, 'w') as f:
        yaml.dump(data, f, sort_keys=False, default_flow_style=False)

    return yaml_path


def gene_family(master_sheet: str, output_dir: str, force: bool = False) -> list:
    output_dir = output_dir or str(Path.cwd())
    logger.info(f"Reading GlobDB protein master sheet: {master_sheet}")

    written = []
    with open(master_sheet, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            gdp_id = row['GDP ID']
            logger.info(f"Building gene family files for {gdp_id}")

            written.append(build_info_yaml(row, output_dir, force=force))

            synteny_path = build_synteny_yaml(row, output_dir, force=force)
            if synteny_path:
                written.append(synteny_path)

            validation_path = build_validation_yaml(row, output_dir, force=force)
            if validation_path:
                written.append(validation_path)

    logger.info(f"Wrote {len(written)} YAML files to {output_dir}")

    return written
