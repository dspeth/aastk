import logging
import shutil
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

from tqdm import tqdm

from .setup import  setup_database
from .schema import TAXONOMY_COLUMNS, HIGH_LEVEL_ENV_COLUMNS, LOW_LEVEL_ENV_COLUMNS
from .sources.sqlite_shards import *
from .sources.tmhmm import *
from .sources.taxonomy import parse_taxonomy_file
from .sources.culture_collection import parse_culture_collection_file
from .sources.environment import parse_environment_file

logger = logging.getLogger(__name__)

def _insert_environment_data(conn, env_path: str, columns: list, label: str):
    logger.info(f"Processing {label} environment data...")
    env_data = parse_environment_file(env_path, columns)
    if env_data:
        column_names = ', '.join(['genome_ID'] + columns)
        placeholders = ', '.join(['?' for _ in range(len(columns) + 1)])
        update_clause = ', '.join([f'{col} = excluded.{col}' for col in columns])
        conn.executemany(f"""
            INSERT INTO genome_data ({column_names})
            VALUES ({placeholders})
            ON CONFLICT(genome_ID) DO UPDATE SET
                {update_clause}
        """, env_data)
        conn.commit()
        logger.info(f"Inserted {label} environment data for {len(env_data)} genomes")

def database(chunk_dir: str = None,
             tmhmm_tar_path: str = None,
             taxonomy_path: str = None,
             culture_collection_path: str = None,
             high_level_environment_path: str = None,
             low_level_environment_path: str = None,
             db_version: str = None,
             output: str = None,
             anvio: bool = False,
             tmhmm: bool = False,
             taxonomy: bool = False,
             culture_collection: bool = False,
             environmental_data: bool = False,
             all_sources: bool = False,
             threads: int = 1,
             force: bool = False):
    if all_sources:
        missing = [
            name for name, path in [
                ('--chunk_dir', chunk_dir),
                ('--tmhmm_tar_path', tmhmm_tar_path),
                ('--taxonomy_path', taxonomy_path),
                ('--culture_collection_path', culture_collection_path),
                ('--high_level_environment_path', high_level_environment_path),
                ('--low_level_environment_path', low_level_environment_path),
            ] if not path
        ]
        if missing:
            raise ValueError(f'--all_sources is True. Missing required path(s): {", ".join(missing)}')

        anvio = True
        tmhmm = True
        taxonomy = True
        culture_collection = True
        environmental_data = True

    output_dir = Path(output) if output else Path('.')
    output_dir.mkdir(parents=True, exist_ok=True)
    db_path = str(output_dir / f"globdb_r{db_version}_aastk.db")

    if force and Path(db_path).exists():
        logger.warning(f"Removing existing database: {db_path}")
        Path(db_path).unlink()

    conn = setup_database(db_path)
    if anvio:
        if not chunk_dir:
            raise ValueError('--anvio is True. --chunk_dir is required.')
        # replace in INSERT with generic call to PROTEIN_COLUMNS, replace VALUES with generic multiplication call
        else:
            BATCH_SIZE = 100
            shard_paths = list(iter_shard_paths(chunk_dir))
            batches = [shard_paths[i:i + BATCH_SIZE] for i in range(0, len(shard_paths), BATCH_SIZE)]

            with ThreadPoolExecutor(max_workers=threads) as executor:
                for batch in tqdm(batches, desc="Processing anvi\'o shards..."):
                    futures = [executor.submit(read_shard, shard_path) for shard_path in batch]
                    records = []
                    for future in as_completed(futures):
                        # read_shard outputs a list of tuples; .extend() flattens these into a new list of tuples
                        records.extend(future.result())
                    conn.executemany("""
                        INSERT INTO protein_data (seqID, parent_ID, aa_length, strand, COG_ID, KEGG_ID, Pfam_ID, cugo_number, protein_seq) 
                        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                        ON CONFLICT(seqID) DO UPDATE SET
                            protein_seq = excluded.protein_seq
                    """, records)
                    conn.commit()

            logger.info("Creating parent_ID index...")
            conn.execute('CREATE INDEX IF NOT EXISTS idx_parent_id ON protein_data(parent_ID)')
            conn.commit()

    if tmhmm is True:
        if not tmhmm_tar_path:
            raise ValueError('--tmhmm is True. -- tmhmm_tar_path must be provided.')
        else:
            logger.info("STEP 3: Processing TMHMM files...")
            tempdir = Path(db_path).parent / 'tmhmm_tempdir'
            if tempdir.exists():
                shutil.rmtree(tempdir)
            tempdir.mkdir(parents=True, exist_ok=True)

            try:
                subprocess.run(["tar", "-xzf", tmhmm_tar_path, "-C", str(tempdir)], check=True)

                tmhmm_files = list(tempdir.rglob("*_tmhmm_clean"))
                logger.info(f"Found {len(tmhmm_files)} TMHMM files")

                for tmhmm_file in tqdm(tmhmm_files, desc="Processing TMHMM"):
                    tmhmm_data = process_tmhmm_file(str(tmhmm_file))
                    if tmhmm_data:
                        conn.executemany("""
                                         INSERT INTO protein_data (seqID, no_tmh)
                                         VALUES (?, ?) ON CONFLICT(seqID) DO
                                         UPDATE SET
                                             no_tmh = excluded.no_tmh
                                         """, [(seqid, no_tmh) for no_tmh, seqid in tmhmm_data])
                        conn.commit()
            finally:
                shutil.rmtree(tempdir)

    if taxonomy:
        if not taxonomy_path:
            raise ValueError('--taxonomy is True. --taxonomy_path must be provided.')
        else:
            logger.info("Processing taxonomy data...")
            taxonomy_data = parse_taxonomy_file(taxonomy_path)
            if taxonomy_data:
                columns = ', '.join(TAXONOMY_COLUMNS)
                placeholders = ', '.join(['?' for _ in TAXONOMY_COLUMNS])
                update_clause = ', '.join(
                    f'{col} = excluded.{col}' for col in TAXONOMY_COLUMNS if col != 'genome_ID'
                )
                conn.executemany(f"""
                    INSERT INTO genome_data ({columns})
                    VALUES ({placeholders})
                    ON CONFLICT(genome_ID) DO UPDATE SET
                        {update_clause}
                """, taxonomy_data)
                conn.commit()
                logger.info(f"Inserted taxonomy data for {len(taxonomy_data)} genomes")

    if culture_collection:
        if not culture_collection_path:
            raise ValueError('--culture_collection is True. --culture_collection_path must be provided.')
        else:
            logger.info("Processing culture collection data...")
            cc_data = parse_culture_collection_file(culture_collection_path)
            if cc_data:
                conn.executemany("""
                    INSERT INTO genome_data (genome_ID, culture_collection)
                    VALUES (?, ?)
                    ON CONFLICT(genome_ID) DO UPDATE SET
                        culture_collection = excluded.culture_collection
                """, cc_data)
                conn.commit()
                logger.info(f"Inserted culture collection data for {len(cc_data)} genomes")

    if environmental_data:
        if not high_level_environment_path or not low_level_environment_path:
            raise ValueError(
                '--environmental_data is True. --high_level_environment_path and '
                '--low_level_environment_path must both be provided.'
            )
        else:
            _insert_environment_data(conn, high_level_environment_path, HIGH_LEVEL_ENV_COLUMNS, 'high-level')
            _insert_environment_data(conn, low_level_environment_path, LOW_LEVEL_ENV_COLUMNS, 'low-level')

    conn.close()

