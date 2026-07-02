import gzip
import sqlite3
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import defaultdict

from ...util import compress_sequence

import logging
logger = logging.getLogger(__name__)

def iter_shard_paths(chunk_dir):
    yield from Path(chunk_dir).rglob("*.db.gz")

def open_shard(shard_path: Path) -> sqlite3.Connection:
    raw = gzip.decompress(shard_path.read_bytes())
    conn = sqlite3.connect(":memory:")
    conn.deserialize(raw)
    return conn

def read_sequences(conn: sqlite3.Connection, prefix: str):
    rows = conn.execute("""
        SELECT gene_callers_id, sequence
        FROM gene_amino_acid_sequences
    """).fetchall()

    return [
        (gene_caller_id, f"{prefix}___{gene_caller_id}", compress_sequence(sequence))
        for gene_caller_id, sequence in rows
    ]

def read_annotations(conn: sqlite3.Connection):
    rows = conn.execute("""
           SELECT gene_callers_id, source, accession, e_value
           FROM gene_functions
           WHERE source IN ('COG24_FUNCTION', 'KOfam', 'Pfam')
        """).fetchall()

    annotations = defaultdict(dict)
    best_pfam_evalue = {}

    for gene_caller_id, source, accession, e_value in rows:
        if source == "COG24_FUNCTION":
            annotations[gene_caller_id]['COG_ID'] = accession
        elif source == "KOfam":
            annotations[gene_caller_id]['KEGG_ID'] = accession
        elif source == "Pfam" and e_value < best_pfam_evalue.get(gene_caller_id, float('inf')):
            best_pfam_evalue[gene_caller_id] = e_value
            annotations[gene_caller_id]['Pfam_ID'] = accession

    return annotations

def read_shard(shard_path: Path):
    prefix = shard_path.name.replace('.db.gz', '')
    conn = open_shard(shard_path)
    try:
        sequences = read_sequences(conn, prefix)
        annotations = read_annotations(conn)

        records = []
        for gene_caller_id, seqID, protein_seq in sequences:
            ann = annotations.get(gene_caller_id, {})
            records.append((
                seqID, ann.get('COG_ID'), ann.get('KEGG_ID'), ann.get('Pfam_ID'), protein_seq
            ))
        return records
    finally:
        conn.close()


