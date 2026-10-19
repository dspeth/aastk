import logging
import sqlite3
from concurrent.futures import ThreadPoolExecutor, as_completed

import matplotlib.pyplot as plt
from tqdm import tqdm
import hashlib
import numpy as np

from aastk.pasr import pasr_select
from aastk.util import ensure_path, determine_dataset_name, read_fasta_to_dict
from aastk.db.schema import (
    BASE_COLUMNS,
    ANNOTATION_COLUMNS,
    TAXONOMY_COLUMNS,
    CULTURE_COLLECTION_COLUMNS,
    HIGH_LEVEL_ENV_COLUMNS,
    LOW_LEVEL_ENV_COLUMNS,
    HIGH_LEVEL_ENV_VIEW,
    LOW_LEVEL_ENV_VIEW,
    HIGH_LEVEL_ENV_CATEGORY_COLUMN,
    LOW_LEVEL_ENV_CATEGORY_COLUMN,
)

logger = logging.getLogger(__name__)

# color maps used for plotting in meta taxonomy plot
COLOR_MAPS = ['Blues', 'Oranges', 'Greens', 'Purples', 'Reds']
DOMAIN_COLORS = {
    'Bacteria': '#3b3b58',
    'Archaea':  '#c9a227',
}

GENOME_ID_MATCH = """
    CASE
        WHEN instr(p.seqID, '___') > 0
        THEN substr(p.seqID, 1, instr(p.seqID, '___') - 1)
        ELSE p.seqID
    END"""


def build_query(include_taxonomy: bool = False,
                include_annotation: bool = False,
                include_culture_collection: bool = False,
                include_high_level_environment: bool = False,
                include_low_level_environment: bool = False,
                include_high_level_env_category: bool = False,
                include_low_level_env_category: bool = False,
                batch_size: int = 1) -> str:
    select_cols = [f'p.{col}' for col in BASE_COLUMNS]

    if include_annotation:
        select_cols.extend([f'p.{col}' for col in ANNOTATION_COLUMNS])

    if include_taxonomy:
        select_cols.extend([f'g.{col}' for col in TAXONOMY_COLUMNS])

    if include_culture_collection:
        select_cols.extend([f'g.{col}' for col in CULTURE_COLLECTION_COLUMNS])

    if include_high_level_environment:
        select_cols.extend([f'g.{col}' for col in HIGH_LEVEL_ENV_COLUMNS])

    if include_low_level_environment:
        select_cols.extend([f'g.{col}' for col in LOW_LEVEL_ENV_COLUMNS])

    if include_high_level_env_category:
        select_cols.append(f'hlenv.{HIGH_LEVEL_ENV_CATEGORY_COLUMN}')

    if include_low_level_env_category:
        select_cols.append(f'llenv.{LOW_LEVEL_ENV_CATEGORY_COLUMN}')

    query = f"SELECT {', '.join(select_cols)} FROM protein_data p"

    needs_genome_join = (
            include_taxonomy or include_culture_collection or
            include_high_level_environment or include_low_level_environment
    )

    if needs_genome_join:
        query += f"""
            LEFT JOIN genome_data g ON
            {GENOME_ID_MATCH} = g.genome_ID
            """

    if include_high_level_env_category:
        query += f"""
            LEFT JOIN {HIGH_LEVEL_ENV_VIEW} hlenv ON
            {GENOME_ID_MATCH} = hlenv.genome_ID
            """

    if include_low_level_env_category:
        query += f"""
            LEFT JOIN {LOW_LEVEL_ENV_VIEW} llenv ON
            {GENOME_ID_MATCH} = llenv.genome_ID
            """

    if batch_size > 1:
        placeholders = ','.join(['?' for _ in range(batch_size)])
        query += f" WHERE p.seqID IN ({placeholders})"
    else:
        query += " WHERE p.seqID = ?"

    return query


QUERY_COLUMNS = {
    **{col: 'p' for col in ANNOTATION_COLUMNS},
    **{col: 'g' for col in TAXONOMY_COLUMNS},
    **{col: 'g' for col in CULTURE_COLLECTION_COLUMNS},
    HIGH_LEVEL_ENV_CATEGORY_COLUMN: 'hlenv',
    LOW_LEVEL_ENV_CATEGORY_COLUMN: 'llenv',
}


def parse_query(query: list) -> dict:
    parsed = {}
    for item in query:
        if '=' not in item:
            raise ValueError(f"Invalid query condition '{item}'. Expected format: column=value")

        column, value = item.split('=', 1)
        column = column.strip()
        value = value.strip()

        if column not in QUERY_COLUMNS:
            raise ValueError(
                f"Unknown metadata column '{column}'. Run 'aastk list_metadata' to view available options."
            )

        parsed.setdefault(column, []).append(value)

    return parsed


def build_metadata_query(query: dict, seq_id_batch: list = None) -> tuple[str, list]:
    tables = {QUERY_COLUMNS[column] for column in query}

    needs_genome_join = 'g' in tables
    needs_hlenv_join = 'hlenv' in tables
    needs_llenv_join = 'llenv' in tables

    sql = "SELECT p.seqID, p.protein_seq FROM protein_data p"

    if needs_genome_join:
        sql += f"""
            LEFT JOIN genome_data g ON
            {GENOME_ID_MATCH} = g.genome_ID
            """

    if needs_hlenv_join:
        sql += f"""
            LEFT JOIN {HIGH_LEVEL_ENV_VIEW} hlenv ON
            {GENOME_ID_MATCH} = hlenv.genome_ID
            """

    if needs_llenv_join:
        sql += f"""
            LEFT JOIN {LOW_LEVEL_ENV_VIEW} llenv ON
            {GENOME_ID_MATCH} = llenv.genome_ID
            """

    conditions = []
    params = []

    for column, values in query.items():
        table = QUERY_COLUMNS[column]
        placeholders = ','.join('?' * len(values))
        conditions.append(f"{table}.{column} IN ({placeholders})")
        params.extend(values)

    if seq_id_batch:
        placeholders = ','.join('?' * len(seq_id_batch))
        conditions.append(f"p.seqID IN ({placeholders})")
        params.extend(seq_id_batch)

    sql += " WHERE " + " AND ".join(conditions)

    return sql, params


def get_header(include_taxonomy: bool,
               include_annotation: bool,
               include_culture_collection: bool,
               include_high_level_environment: bool,
               include_low_level_environment: bool,
               include_high_level_env_category: bool = False,
               include_low_level_env_category: bool = False):

    header = BASE_COLUMNS.copy()

    if include_annotation:
        header.extend(ANNOTATION_COLUMNS)

    if include_taxonomy:
        header.extend(TAXONOMY_COLUMNS)

    if include_culture_collection:
        header.extend(CULTURE_COLLECTION_COLUMNS)

    if include_high_level_environment:
        header.extend(HIGH_LEVEL_ENV_COLUMNS)

    if include_low_level_environment:
        header.extend(LOW_LEVEL_ENV_COLUMNS)

    if include_high_level_env_category:
        header.append(HIGH_LEVEL_ENV_CATEGORY_COLUMN)

    if include_low_level_env_category:
        header.append(LOW_LEVEL_ENV_CATEGORY_COLUMN)

    return header


# one representative column per category: presence of data there indicates
# that the corresponding ingestion step (build.py) has actually been run
AVAILABILITY_CHECKS = {
    'annotation': ('protein_data', 'COG_ID'),
    'taxonomy': ('genome_data', 'domain'),
    'culture_collection': ('genome_data', 'culture_collection'),
    'high_level_environment': ('genome_data', 'animal_associated'),
    'low_level_environment': ('genome_data', 'human'),
    'high_level_env_category': (HIGH_LEVEL_ENV_VIEW, HIGH_LEVEL_ENV_CATEGORY_COLUMN),
    'low_level_env_category': (LOW_LEVEL_ENV_VIEW, LOW_LEVEL_ENV_CATEGORY_COLUMN),
}


def check_metadata_availability(db_path: str,
                                include_taxonomy: bool,
                                include_annotation: bool,
                                include_culture_collection: bool,
                                include_high_level_environment: bool,
                                include_low_level_environment: bool,
                                include_high_level_env_category: bool = False,
                                include_low_level_env_category: bool = False) -> list:
    requested = {
        'annotation': include_annotation,
        'taxonomy': include_taxonomy,
        'culture_collection': include_culture_collection,
        'high_level_environment': include_high_level_environment,
        'low_level_environment': include_low_level_environment,
        'high_level_env_category': include_high_level_env_category,
        'low_level_env_category': include_low_level_env_category,
    }

    unavailable = []
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    for category, is_requested in requested.items():
        if not is_requested:
            continue
        table, column = AVAILABILITY_CHECKS[category]
        cursor.execute(f"SELECT 1 FROM {table} WHERE {column} IS NOT NULL LIMIT 1")
        if cursor.fetchone() is None:
            unavailable.append(category)

    conn.close()
    return unavailable


def process_batch(db_path: str,
                  batch: list,
                  include_taxonomy: bool,
                  include_annotation: bool,
                  include_culture_collection: bool,
                  include_high_level_environment: bool,
                  include_low_level_environment: bool,
                  include_high_level_env_category: bool = False,
                  include_low_level_env_category: bool = False):

    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    query = build_query(
        include_taxonomy=include_taxonomy,
        include_annotation=include_annotation,
        include_culture_collection=include_culture_collection,
        include_high_level_environment=include_high_level_environment,
        include_low_level_environment=include_low_level_environment,
        include_high_level_env_category=include_high_level_env_category,
        include_low_level_env_category=include_low_level_env_category,
        batch_size=len(batch)
    )

    cursor.execute(query, batch)
    results = cursor.fetchall()

    result_dict = {row[0]: row for row in results}

    conn.close()
    return result_dict


def meta(db_path: str,
         fasta: str,
         id_list: str,
         output: str,
         threads: int = 1,
         include_annotation: bool = False,
         include_taxonomy: bool = False,
         include_culture_collection: bool = False,
         include_high_level_environment: bool = False,
         include_low_level_environment: bool = False,
         include_high_level_env_category: bool = False,
         include_low_level_env_category: bool = False,
         high_level_env_stacked_bar_plot: bool = False,
         low_level_env_stacked_bar_plot: bool = False,
         taxonomy_stacked_bar_plot: bool = False,
         all_metadata: bool = False,
         force: bool = False):

    if fasta:
        seq_dict = read_fasta_to_dict(fasta)
        seq_ids = list(seq_dict.keys())
        dataset_name = determine_dataset_name(fasta, '.', 0)
        output_path = ensure_path(
            output,
            f'{dataset_name}_metadata.tsv',
            force=force
        )
    elif id_list:
        with open(id_list, 'r') as f:
            seq_ids = [line.strip() for line in f]
        dataset_name = determine_dataset_name(id_list, '.', 0)
        output_path = ensure_path(
            output,
            f'{dataset_name}_metadata.tsv',
            force=force
        )
    else:
        logger.error(
            "No valid input file found. Please specify path to protein FASTA file OR sequence ID list (one ID per line)."
        )
        return

    if include_high_level_env_category is True:
        high_level_env_stacked_bar_plot = True
    elif include_low_level_env_category is True:
        low_level_env_stacked_bar_plot = True
    elif include_taxonomy is True:
        taxonomy_stacked_bar_plot = True


    if all_metadata:
        include_annotation = True
        include_taxonomy = True
        include_culture_collection = True
        include_high_level_environment = True
        include_low_level_environment = True
        include_high_level_env_category = True
        include_low_level_env_category = True
        high_level_env_stacked_bar_plot = True
        low_level_env_stacked_bar_plot = True
        taxonomy_stacked_bar_plot = True

    if not any([
        include_annotation,
        include_taxonomy,
        include_culture_collection,
        include_high_level_environment,
        include_low_level_environment,
        include_high_level_env_category,
        include_low_level_env_category
    ]):
        logger.error(
            "Must specify at least one metadata type "
            "(--taxonomy, --annotation, --culture-collection, "
            "--high-level-environment, --low-level-environment, "
            "--high-level-env-category, --low-level-env-category, or --all-metadata)"
        )
        return

    unavailable = check_metadata_availability(
        db_path,
        include_taxonomy,
        include_annotation,
        include_culture_collection,
        include_high_level_environment,
        include_low_level_environment,
        include_high_level_env_category,
        include_low_level_env_category
    )
    if unavailable:
        logger.error(
            f"Requested metadata not available in database: {', '.join(unavailable)}"
        )
        return

    BATCH_SIZE = 500

    header = get_header(
        include_taxonomy,
        include_annotation,
        include_culture_collection,
        include_high_level_environment,
        include_low_level_environment,
        include_high_level_env_category,
        include_low_level_env_category
    )

    # Split seq_ids into batches
    batches = [
        seq_ids[i:i + BATCH_SIZE]
        for i in range(0, len(seq_ids), BATCH_SIZE)
    ]

    missing_ids = set()
    results = {}

    # Process batches in parallel
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {
            executor.submit(
                process_batch,
                db_path,
                batch,
                include_taxonomy,
                include_annotation,
                include_culture_collection,
                include_high_level_environment,
                include_low_level_environment,
                include_high_level_env_category,
                include_low_level_env_category
            ): idx
            for idx, batch in enumerate(batches)
        }

        for future in tqdm(
                as_completed(futures),
                total=len(batches),
                desc="Querying database"
        ):
            try:
                results.update(future.result())
            except Exception as e:
                logger.error(
                    f"Error processing batch {futures[future]}: {e}"
                )

    # Write results in original order
    with open(output_path, 'w') as out:
        out.write('\t'.join(header) + '\n')

        for seq_id in seq_ids:
            if seq_id in results:
                row = [
                    str(x) if x is not None else ''
                    for x in results[seq_id]
                ]
                out.write('\t'.join(row) + '\n')
            else:
                missing_ids.add(seq_id)

    logger.info(f"Metadata has been saved to {output_path}")

    if include_taxonomy is True:
        plot_taxonomy(header, results, output, dataset_name, force=force)

    if missing_ids:
        logger.warning(
            f"{len(missing_ids)} sequence IDs were not found in the database"
        )

def taxon_color(header: list,
                results: dict,
                top_n_dict: dict,
                top_n: int):
    phylum_to_colormap = {}
    for i, phylum in enumerate(top_n_dict['phylum'].keys()):
        if phylum != 'Other':
            phylum_to_colormap[phylum] = COLOR_MAPS[i]

    # save top_n shades for each color map; 1 to 0.25 to begin with strongest shade
    color_shades = {}
    for color_map in COLOR_MAPS:
        cmap = plt.get_cmap(color_map)
        # get RGB array for equidistantly spaced points on the Colormap scale
        color_shades[color_map] = cmap(np.linspace(1, 0.25, top_n + 1))

    # map phylum to color_shade
    phylum_color = {phylum: color_shades[phylum_to_colormap[phylum]][0] for phylum in phylum_to_colormap.keys()}

    # map taxon to phylum and assign color shade according to top_n directed sort
    phylum_index = list(header).index('phylum')
    taxon_to_phylum = {}
    # only map phylum to ranks below phylum
    for rank, counts in list(top_n_dict.items())[2:]:
        rank_index = header.index(rank)
        taxon_to_phylum[rank] = {}
        for row in results.values():
            taxon = row[rank_index]
            if taxon in counts:
                taxon_to_phylum[rank][taxon] = row[phylum_index]

    # map taxon to correct phylum-derived shade
    taxon_to_phylum_color_shade = {}
    for rank, counts in list(top_n_dict.items())[2:]:
        taxon_to_phylum_color_shade[rank] = {}
        phylum_count = {}
        for taxon in counts:
            if taxon == 'Other':
                continue
            phylum = taxon_to_phylum[rank][taxon]
            # keep track of how many times phylum was encountered starting with 1
            # as 0 is occupied by phylum itself when it comes to color mapping
            phylum_count[phylum] = phylum_count.get(phylum, 0) + 1
            shades = color_shades[phylum_to_colormap[phylum]]
            taxon_to_phylum_color_shade[rank][taxon] = shades[phylum_count[phylum]]

    return phylum_color, taxon_to_phylum_color_shade


def plot_taxonomy(header: list,
                  results: dict,
                  output: str,
                  dataset_name: str,
                  top_n: int = 5,
                  force: bool = False):
    # start with phylum; omit genome_ID, genus and species
    taxonomy_dict = {key: {} for key in TAXONOMY_COLUMNS[1:6]}
    for column in taxonomy_dict.keys():
        index = next((id for id, item in enumerate(header) if column in item), None)
        for key, value in results.items():
            if isinstance(value, (tuple, list)) and index < len(value):
                if results[key][index] not in taxonomy_dict[column].keys():
                    taxonomy_dict[column][results[key][index]] = 1
                else:
                    taxonomy_dict[column][results[key][index]] += 1
            else:
                logger.warning(f"Key '{key}' has incomplete data: {value}")

    # get top n per taxonomic rank
    # outer loop: iterate over ranks -> column = rank name, counts = {taxon: count}
    # inner sort: x = (taxon, count) tuple, x[1] = count -> sort taxa by count, descending
    # slice the sorted list before converting back to dict (dicts can't be sliced)
    top_n_dict = {
    column: dict(sorted(counts.items(), key=lambda x: x[1], reverse=True)[:top_n])
    for column, counts in taxonomy_dict.items()
    }

    total_sequences = len(results)
    for key in top_n_dict.keys():
        top_n_sum = sum(top_n_dict[key].values())
        other_count = total_sequences - top_n_sum
        if other_count > 0:
            top_n_dict[key]['Other'] = other_count

    phylum_color, taxon_shades = taxon_color(header, results, top_n_dict, top_n)

    fig, axes = plt.subplots(1, len(top_n_dict), figsize=(3.5 * len(top_n_dict), 8), sharey=True)
    for ax, (column, counts) in zip(axes, top_n_dict.items()):
        rank_name = 'Order' if column == 'order_tax' else column.capitalize()
        bottom = 0
        for taxon, count in counts.items():
            if column == 'domain':
                color = DOMAIN_COLORS.get(str(taxon).replace('d__', ''), '#cccccc')
            elif column == 'phylum':
                color = phylum_color.get(taxon, '#cccccc')
            else:
                color = taxon_shades.get(column, {}).get(taxon, '#cccccc')
            ax.bar(0, count, bottom=bottom, width=0.8, label=str(taxon), color=color, linewidth=0)
            bottom += count
        ax.set_xlim(-0.75, 0.75)
        ax.set_ylim(0, total_sequences)
        ax.set_xticks([])
        ax.set_title(rank_name)
        for side in ['top', 'right', 'bottom']:
            ax.spines[side].set_visible(False)
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[::-1], labels[::-1], loc='upper center', bbox_to_anchor=(0.5, -0.02), frameon=False, fontsize=8)
    axes[0].set_ylabel('Number of sequences')

    plot_path = ensure_path(
        output,
        f'{dataset_name}_taxonomy_top{top_n}_stacked.svg',
        force=force
    )
    plt.tight_layout()
    plt.savefig(plot_path, format='svg', bbox_inches='tight')
    plt.close(fig)
    logger.info(f"Plot saved to {plot_path}")

def meta_plot(results: dict,
              high_level_env_stacked_bar_plot: bool = False,
              low_level_env_stacked_bar_plot: bool = False,
              taxonomy_stacked_bar_plot: bool = False
              ):
    plot_taxonomy(results)


def list_metadata():
    print("───────────────────────────")
    print("Protein metadata")
    print("───────────────────────────")
    print("  • sequence information")
    for item in BASE_COLUMNS[1:]:
        print(f"      - {item}")
    print()
    print("  • annotations")
    for item in ANNOTATION_COLUMNS:
        print(f"      - {item}")
    print()  # blank line

    print("───────────────────────────")
    print("Genome metadata")
    print("───────────────────────────")
    print("  • taxonomy")
    for item in TAXONOMY_COLUMNS:
        print(f"      - {item}")
    print()

    print("  • culture information")
    for item in CULTURE_COLLECTION_COLUMNS:
        print(f"      - {item}")
    print()

    print("  • high-level environment data")
    for item in HIGH_LEVEL_ENV_COLUMNS:
        print(f"      - {item}")
    print()

    print("  • low-level environment data")
    for item in LOW_LEVEL_ENV_COLUMNS:
        print(f"      - {item}")
    print()

    print("  • derived environment category (dominant category, >=50%, else 'diverse')")
    print(f"      - {HIGH_LEVEL_ENV_CATEGORY_COLUMN}")
    print(f"      - {LOW_LEVEL_ENV_CATEGORY_COLUMN}")
    print()
