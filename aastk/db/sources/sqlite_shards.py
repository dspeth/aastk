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
        SELECT s.gene_callers_id, s.sequence
        FROM gene_amino_acid_sequences s
        JOIN genes_in_contigs g ON g.gene_callers_id = s.gene_callers_id
        WHERE g.call_type = 1
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



def read_contig_info(conn: sqlite3.Connection):
    rows = conn.execute("""
        SELECT gene_callers_id, contig, start, stop, direction
        FROM genes_in_contigs
        WHERE call_type = 1
        """).fetchall()

    contig_info = defaultdict(dict)

    for gene_caller_id, contig, start, stop, direction in rows:
        contig_info[gene_caller_id]['contig'] = contig
        contig_info[gene_caller_id]['start'] = start
        contig_info[gene_caller_id]['stop'] = stop
        contig_info[gene_caller_id]['aa_length'] = int((abs(int(stop) - int(start)) + 1)/3)
        if direction == 'f':
            contig_info[gene_caller_id]['direction'] = '+'
        elif direction == 'r':
            contig_info[gene_caller_id]['direction'] = '-'

    sorted_contig_info = {k: v for k, v in sorted(contig_info.items(),
                          key=lambda x: (x[1]['contig'], x[1]['start']))}

    cugo_number = 0
    previous_key = None
    for key in sorted_contig_info.keys():
        if previous_key is None or (sorted_contig_info[key]['direction'] == sorted_contig_info[previous_key]['direction'] and sorted_contig_info[key]['contig'] == sorted_contig_info[previous_key]['contig']):
            contig_info[key]['cugo_number'] = cugo_number
            previous_key = key
        else:
            cugo_number += 1
            contig_info[key]['cugo_number'] = cugo_number
            previous_key = key

    return contig_info

def read_shard(shard_path: Path):
    prefix = shard_path.name.replace('.db.gz', '')
    conn = open_shard(shard_path)
    try:
        sequences = read_sequences(conn, prefix)
        annotations = read_annotations(conn)
        contig_info = read_contig_info(conn)

        records = []
        for gene_caller_id, seqID, protein_seq in sequences:
            ann = annotations.get(gene_caller_id, {})
            contig = contig_info.get(gene_caller_id, {})

            records.append((
                seqID, contig.get('contig'), contig.get('aa_length'), contig.get('direction'), ann.get('COG_ID', 'NA'),
                ann.get('KEGG_ID', 'NA'), ann.get('Pfam_ID', 'NA'), contig.get('cugo_number'), protein_seq
            ))
        return records
    finally:
        conn.close()


