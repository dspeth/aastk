import sqlite3
from concurrent.futures import ProcessPoolExecutor, as_completed

from aastk.util import ensure_path, decompress_sequence


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
