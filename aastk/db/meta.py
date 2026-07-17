import logging
import sqlite3
from concurrent.futures import ThreadPoolExecutor, as_completed

from tqdm import tqdm

from aastk.util import ensure_path, determine_dataset_name, read_fasta_to_dict
from aastk.db.schema import (
    BASE_COLUMNS,
    ANNOTATION_COLUMNS,
    TAXONOMY_COLUMNS,
    CULTURE_COLLECTION_COLUMNS,
    HIGH_LEVEL_ENV_COLUMNS,
    LOW_LEVEL_ENV_COLUMNS,
)

logger = logging.getLogger(__name__)


def build_query(include_taxonomy: bool = False,
                include_annotation: bool = False,
                include_culture_collection: bool = False,
                include_high_level_environment: bool = False,
                include_low_level_environment: bool = False,
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

    query = f"SELECT {', '.join(select_cols)} FROM protein_data p"

    needs_genome_join = (
            include_taxonomy or include_culture_collection or
            include_high_level_environment or include_low_level_environment
    )

    if needs_genome_join:
        query += """
            LEFT JOIN genome_data g ON
            CASE
                WHEN instr(p.seqID, '___') > 0
                THEN substr(p.seqID, 1, instr(p.seqID, '___') - 1)
                ELSE p.seqID
            END = g.genome_ID
            """

    if batch_size > 1:
        placeholders = ','.join(['?' for _ in range(batch_size)])
        query += f" WHERE p.seqID IN ({placeholders})"
    else:
        query += " WHERE p.seqID = ?"

    return query


def get_header(include_taxonomy: bool,
               include_annotation: bool,
               include_culture_collection: bool,
               include_high_level_environment: bool,
               include_low_level_environment: bool):

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

    return header


# one representative column per category: presence of data there indicates
# that the corresponding ingestion step (build.py) has actually been run
AVAILABILITY_CHECKS = {
    'annotation': ('protein_data', 'COG_ID'),
    'taxonomy': ('genome_data', 'domain'),
    'culture_collection': ('genome_data', 'culture_collection'),
    'high_level_environment': ('genome_data', 'animal_associated'),
    'low_level_environment': ('genome_data', 'human'),
}


def check_metadata_availability(db_path: str,
                                include_taxonomy: bool,
                                include_annotation: bool,
                                include_culture_collection: bool,
                                include_high_level_environment: bool,
                                include_low_level_environment: bool) -> list:
    requested = {
        'annotation': include_annotation,
        'taxonomy': include_taxonomy,
        'culture_collection': include_culture_collection,
        'high_level_environment': include_high_level_environment,
        'low_level_environment': include_low_level_environment,
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
                  include_low_level_environment: bool):

    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    query = build_query(
        include_taxonomy=include_taxonomy,
        include_annotation=include_annotation,
        include_culture_collection=include_culture_collection,
        include_high_level_environment=include_high_level_environment,
        include_low_level_environment=include_low_level_environment,
        batch_size=len(batch)
    )

    cursor.execute(query, batch)
    results = cursor.fetchall()

    result_dict = {row[0]: row for row in results}

    conn.close()
    return batch, result_dict


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

    # Enable everything if --all-metadata is set
    if all_metadata:
        include_annotation = True
        include_taxonomy = True
        include_culture_collection = True
        include_high_level_environment = True
        include_low_level_environment = True

    # Ensure at least one metadata type is requested
    if not any([
        include_annotation,
        include_taxonomy,
        include_culture_collection,
        include_high_level_environment,
        include_low_level_environment
    ]):
        logger.error(
            "Must specify at least one metadata type "
            "(--taxonomy, --annotation, --culture-collection, "
            "--high-level-environment, --low-level-environment, or --all-metadata)"
        )
        return

    unavailable = check_metadata_availability(
        db_path,
        include_taxonomy,
        include_annotation,
        include_culture_collection,
        include_high_level_environment,
        include_low_level_environment
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
        include_low_level_environment
    )

    # Split seq_ids into batches
    batches = [
        seq_ids[i:i + BATCH_SIZE]
        for i in range(0, len(seq_ids), BATCH_SIZE)
    ]

    found_ids = set()
    missing_ids = set()
    batch_results = {}

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
                include_low_level_environment
            ): idx
            for idx, batch in enumerate(batches)
        }

        for future in tqdm(
                as_completed(futures),
                total=len(batches),
                desc="Querying database"
        ):
            batch_idx = futures[future]
            try:
                batch, result_dict = future.result()
                batch_results[batch_idx] = (batch, result_dict)
            except Exception as e:
                logger.error(
                    f"Error processing batch {batch_idx}: {e}"
                )

    # Write results in original order
    with open(output_path, 'w') as out:
        out.write('\t'.join(header) + '\n')

        for batch_idx in sorted(batch_results.keys()):
            batch, result_dict = batch_results[batch_idx]

            for seq_id in batch:
                if seq_id in result_dict:
                    found_ids.add(seq_id)
                    row = [
                        str(x) if x is not None else ''
                        for x in result_dict[seq_id]
                    ]
                    out.write('\t'.join(row) + '\n')
                else:
                    missing_ids.add(seq_id)

    logger.info(f"Metadata has been saved to {output_path}")

    if missing_ids:
        logger.warning(
            f"{len(missing_ids)} sequence IDs were not found in the database"
        )


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
