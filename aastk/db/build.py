import logging
import shutil
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

from tqdm import tqdm

from .setup import  setup_database
from .sources.sqlite_shards import *
from .sources.tmhmm import *

logger = logging.getLogger(__name__)

def database(chunk_dir: str,
             tmhmm_tar_path: str,
             db_path: str,
             anvio: bool = False,
             tmhmm: bool = False,
             taxonomy: bool = False,
             environmental_data: bool = False,
             all_sources: bool = False,
             threads: int = 1):
    conn = setup_database(db_path)

    if all_sources:
        anvio = True
        tmhmm = True
        taxonomy = True
        environmental_data = True
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
                        records.extend(future.result())
                    conn.executemany("""
                        INSERT INTO protein_data (seqID, parent_ID, aa_length, strand, COG_ID, KEGG_ID, Pfam_ID, cugo_number, protein_seq) 
                        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                        ON CONFLICT(seqID) DO UPDATE SET
                            protein_seq = excluded.protein_seq
                    """, records)
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

    conn.close()

