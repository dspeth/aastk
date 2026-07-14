import sqlite3

from aastk.db.schema import PROTEIN_SCHEMA, GENOME_SCHEMA, PROTEIN_COLUMNS


def setup_database(db_path: str) -> sqlite3.Connection:
    conn = sqlite3.connect(db_path)

    protein_cols = ',\n                     '.join(
        [f'{col} {dtype}' for col, dtype in PROTEIN_SCHEMA.items()]
    )

    conn.execute(f'''
                 CREATE TABLE IF NOT EXISTS protein_data (
                    {protein_cols}
                 )
                 ''')

    genome_cols = ',\n                     '.join(
        [f'{col} {dtype}' for col, dtype in GENOME_SCHEMA.items()]
    )

    conn.execute(f'''
                 CREATE TABLE IF NOT EXISTS genome_data (
                    {genome_cols}
                 )
                 ''')

    view_cols = [col for col in PROTEIN_COLUMNS if col != 'protein_seq']
    view_select = ',\n                     '.join(view_cols)

    conn.execute(f'''
                 CREATE VIEW IF NOT EXISTS protein_data_readable AS
                 SELECT
                     {view_select},
                     CASE WHEN protein_seq IS NULL THEN NULL ELSE '<COMPRESSED_BLOB>' END as protein_seq_status
                 FROM protein_data
                 ''')


    conn.commit()
    return conn
