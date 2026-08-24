import json
import logging
import sqlite3
from concurrent.futures import ProcessPoolExecutor, as_completed

from tqdm import tqdm

from aastk.util import ensure_path, decompress_sequence, determine_dataset_name, read_fasta_to_dict
from aastk.db.meta import build_metadata_query

logger = logging.getLogger(__name__)


def stream_sequence_ids(db_path, batch_size):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    cursor.execute("SELECT seqID FROM protein_data")

    batch = []
    for row in cursor:
        batch.append(row[0])
        if len(batch) >= batch_size:
            yield batch
            batch = []
    if batch:
        yield batch

    conn.close()


def fetch_sequences(db_path, seq_ids):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    placeholders = ', '.join(['?'] * len(seq_ids))
    cursor.execute(f"""SELECT seqID, protein_seq FROM protein_data WHERE seqID in ({placeholders})""", seq_ids)

    records = []
    for seq_id, blob in cursor.fetchall():
        if blob:
            records.append((seq_id, decompress_sequence(blob)))

    return records


def export_fasta(db_path: str,
                  output: str,
                  threads: int = 1,
                  force: bool = False):
    version = db_path.split('_')[1].replace('r', '')
    protein_fasta_file = ensure_path(output, f'globdb_r{version}_all_prot.faa', force=force)

    batch_size = 900

    with open(protein_fasta_file, 'w') as out:
        with ProcessPoolExecutor(max_workers=threads) as executor:
            futures = set()

            for batch in stream_sequence_ids(db_path, batch_size):
                futures.add(executor.submit(fetch_sequences, db_path, batch))

                if len(futures) >= threads:
                    done = next(as_completed(futures))
                    futures.remove(done)

                    for seqID, seq in done.result():
                        out.write(f">{seqID}\n{seq}\n")

            for future in as_completed(futures):
                for seqID, seq in future.result():
                    out.write(f">{seqID}\n{seq}\n")

    return protein_fasta_file


def fetch_metadata_batch(db_path: str, query: dict, batch: list):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    sql, params = build_metadata_query(query, seq_id_batch=batch)
    cursor.execute(sql, params)

    records = []
    for seqID, blob in cursor.fetchall():
        if blob:
            records.append((seqID, decompress_sequence(blob)))

    conn.close()
    return records


def export_by_metadata(db_path: str,
                        query: dict,
                        output: str,
                        fasta: str = None,
                        id_list: str = None,
                        threads: int = 1,
                        force: bool = False) -> str:
    dataset_name = determine_dataset_name(db_path, '.', 0)
    output_path = ensure_path(output, f'{dataset_name}_export.faa', force=force)
    query_info_path = ensure_path(output, f'{dataset_name}_export_query.json', force=force)

    written = 0

    if fasta or id_list:
        if fasta:
            seq_ids = list(read_fasta_to_dict(fasta).keys())
        else:
            with open(id_list, 'r') as f:
                seq_ids = [line.strip() for line in f]

        batch_size = 500
        batches = [seq_ids[i:i + batch_size] for i in range(0, len(seq_ids), batch_size)]

        with open(output_path, 'w') as out:
            with ProcessPoolExecutor(max_workers=threads) as executor:
                futures = {
                    executor.submit(fetch_metadata_batch, db_path, query, batch): idx
                    for idx, batch in enumerate(batches)
                }

                for future in tqdm(as_completed(futures), total=len(batches), desc="Querying database"):
                    for seqID, seq in future.result():
                        out.write(f">{seqID}\n{seq}\n")
                        written += 1
    else:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()

        sql, params = build_metadata_query(query)
        cursor.execute(sql, params)

        with open(output_path, 'w') as out:
            for seqID, blob in cursor:
                if blob:
                    out.write(f">{seqID}\n{decompress_sequence(blob)}\n")
                    written += 1

        conn.close()

    with open(query_info_path, 'w') as f:
        json.dump({
            'db_path': db_path,
            'query': query,
            'fasta': fasta,
            'id_list': id_list,
            'sequences_written': written,
        }, f, indent=2)

    logger.info(f"Exported {written} sequences to {output_path}")
    logger.info(f"Query recorded at {query_info_path}")

    return output_path, query_info_path
