import sqlite3

from ..util import ensure_path
from .schema import TAXONOMY_COLUMNS, HIGH_LEVEL_ENV_COLUMNS, LOW_LEVEL_ENV_COLUMNS


def database_check(db_path: str,
                   output: str,
                   force: bool = False):
    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()

        stages = [
            (
                "COG",
                """
                SELECT seqID FROM protein_data
                WHERE COG_ID IS NULL
                """
            ),
            (
                "TMHMM",
                """
                SELECT seqID FROM protein_data
                WHERE no_tmh IS NULL
                """
            ),
            (
                "KEGG",
                """
                SELECT seqID FROM protein_data
                WHERE KEGG_ID IS NULL
                """
            ),
            (
                "Pfam",
                """
                SELECT seqID FROM protein_data
                WHERE Pfam_ID IS NULL
                """
            ),
            (
                "Protein_sequences",
                """
                SELECT seqID FROM protein_data
                WHERE protein_seq IS NULL
                """
            ),
        ]

        for stage, query in stages:
            print(f"Checking for missing data: {stage}")
            cursor.execute(query)
            missing = [row[0] for row in cursor.fetchall()]

            out_file = ensure_path(output, f"missing_{stage}.txt", force=force)
            with open(out_file, 'w') as f:
                for seqid in missing:
                    f.write(f"{seqid}\n")

        genome_stages = [
            (
                "Taxonomy",
                f"""
                SELECT genome_ID FROM genome_data
                WHERE {" OR ".join(f"{col} IS NULL" for col in TAXONOMY_COLUMNS)}
                """
            ),
            (
                "Culture_collection",
                """
                SELECT genome_ID FROM genome_data
                WHERE culture_collection IS NULL
                """
            ),
            (
                "High_level_environment",
                f"""
                SELECT genome_ID FROM genome_data
                WHERE {" OR ".join(f"{col} IS NULL" for col in HIGH_LEVEL_ENV_COLUMNS)}
                """
            ),
            (
                "Low_level_environment",
                f"""
                SELECT genome_ID FROM genome_data
                WHERE {" OR ".join(f"{col} IS NULL" for col in LOW_LEVEL_ENV_COLUMNS)}
                """
            ),
        ]

        for stage, query in genome_stages:
            cursor.execute(query)
            missing = [row[0] for row in cursor.fetchall()]

            out_file = ensure_path(output, f"missing_{stage}.txt", force=force)
            with open(out_file, 'w') as f:
                for genome_id in missing:
                    f.write(f"{genome_id}\n")
