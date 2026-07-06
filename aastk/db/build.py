from .setup import  setup_database
from .sources.sqlite_shards import *

def database(chunk_dir, db_path):
    conn = setup_database(db_path)

    # replace in INSERT with generic call to PROTEIN_COLUMNS, replace VALUES with generic multiplication call
    for shard_path in iter_shard_paths(chunk_dir):
        records = read_shard(shard_path)
        conn.executemany("""
            INSERT INTO protein_data (seqID, parent_ID, aa_length, strand, COG_ID, KEGG_ID, Pfam_ID, protein_seq) 
            VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(seqID) DO UPDATE SET
                protein_seq = excluded.protein_seq
        """, records)
        conn.commit()

    conn.close()

