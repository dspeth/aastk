import sqlite3

from aastk.db.schema import (
    PROTEIN_SCHEMA,
    GENOME_SCHEMA,
    PROTEIN_COLUMNS,
    HIGH_LEVEL_ENV_COLUMNS,
    LOW_LEVEL_ENV_COLUMNS,
    HIGH_LEVEL_ENV_VIEW,
    LOW_LEVEL_ENV_VIEW,
    HIGH_LEVEL_ENV_CATEGORY_COLUMN,
    LOW_LEVEL_ENV_CATEGORY_COLUMN,
)


def _env_category_view_sql(view_name: str, columns: list, category_column: str) -> str:
    col_list = ', '.join(columns)
    tie_expr = ' + '.join(f'({col} = max_val)' for col in columns)
    case_whens = '\n                '.join(f"WHEN {col} = max_val THEN '{col}'" for col in columns)

    return f'''
                 CREATE VIEW IF NOT EXISTS {view_name} AS
                 WITH env AS (
                     SELECT genome_ID,
                            {col_list},
                            MAX({col_list}) AS max_val
                     FROM genome_data
                 )
                 SELECT genome_ID,
                     CASE
                         WHEN max_val IS NULL THEN NULL
                         WHEN max_val < 0.5 THEN 'diverse'
                         WHEN ({tie_expr}) > 1 THEN 'diverse'
                         {case_whens}
                     END AS {category_column}
                 FROM env
                 '''


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

    conn.execute(_env_category_view_sql(HIGH_LEVEL_ENV_VIEW, HIGH_LEVEL_ENV_COLUMNS, HIGH_LEVEL_ENV_CATEGORY_COLUMN))
    conn.execute(_env_category_view_sql(LOW_LEVEL_ENV_VIEW, LOW_LEVEL_ENV_COLUMNS, LOW_LEVEL_ENV_CATEGORY_COLUMN))

    conn.commit()
    return conn
