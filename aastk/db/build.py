from .setup import  setup_database
from .sources.sqlite_shards import *

def database(chunk_dir, db_path):
    conn = setup_database(db_path)

    for shard_path in iter_shard_paths(chunk_dir):
        records = read_shard(shard_path)
        conn.executemany("""
            INSERT INTO protein_data (seqID, COG_ID, KEGG_ID, Pfam_ID, protein_seq)
            VALUES (?, ?, ?, ?, ?)
            ON CONFLICT(seqID) DO UPDATE SET
                protein_seq = excluded.protein_seq
        """, records)
        conn.commit()

    conn.close()

