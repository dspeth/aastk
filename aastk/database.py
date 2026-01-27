from .cugo import *
from .util import *

import sqlite3
import subprocess
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed
from textwrap import dedent


logger = logging.getLogger(__name__)

# define database table schema
PROTEIN_SCHEMA = {
    'seqID': 'TEXT PRIMARY KEY',
    'parent_ID': 'TEXT',
    'aa_length': 'INTEGER',
    'strand': 'TEXT',
    'COG_ID': 'TEXT',
    'KEGG_ID': 'TEXT',
    'Pfam_ID': 'TEXT',
    'cugo_number': 'INTEGER',
    'no_tmh': 'INTEGER',
    'protein_seq': 'BLOB'
}

GENOME_SCHEMA = {
    'genome_ID': 'TEXT PRIMARY KEY',
    'domain': 'TEXT',
    'phylum': 'TEXT',
    'class': 'TEXT',
    'order_tax': 'TEXT',
    'family': 'TEXT',
    'genus': 'TEXT',
    'species': 'TEXT',
    'culture_collection': 'INT',
    'animal_associated': 'REAL',
    'aquatic': 'REAL',
    'built': 'REAL',
    'other': 'REAL',
    'sediment': 'REAL',
    'soil': 'REAL',
    'unassigned_high_level': 'REAL',
    'human': 'REAL',
    'invertebrate': 'REAL',
    'other_vertebrate': 'REAL',
    'unspecified_animal': 'REAL',
    'aquatic_other': 'REAL',
    'aquatic_unspecified': 'REAL',
    'freshwater': 'REAL',
    'groundwater': 'REAL',
    'marine': 'REAL',
    'built_other': 'REAL',
    'drinking_water': 'REAL',
    'wastewater': 'REAL',
    'air': 'REAL',
    'bacteria': 'REAL',
    'eukaryote_other': 'REAL',
    'food': 'REAL',
    'geothermal': 'REAL',
    'hypersaline': 'REAL',
    'other_unspecified': 'REAL',
    'plant_associated': 'REAL',
    'subsurface': 'REAL',
    'synthetic': 'REAL',
    'viral': 'REAL',
    'freshwater_sediment': 'REAL',
    'marine_sediment': 'REAL',
    'sediment_unspecified': 'REAL',
    'desert': 'REAL',
    'forest': 'REAL',
    'rhizosphere': 'REAL',
    'soil_agricultural': 'REAL',
    'soil_other': 'REAL',
    'soil_unspecified': 'REAL',
    'tundra_wetland': 'REAL',
    'unassigned_low_level': 'REAL'
}

# extract column name lists for queries (same order)
PROTEIN_COLUMNS = list(PROTEIN_SCHEMA.keys())
GENOME_COLUMNS = list(GENOME_SCHEMA.keys())

# categorize columns for later use
BASE_COLUMNS = ['seqID', 'parent_ID', 'aa_length', 'strand']

ANNOTATION_COLUMNS = ['COG_ID', 'KEGG_ID', 'Pfam_ID']

TAXONOMY_COLUMNS = [
    'genome_ID', 'domain', 'phylum', 'class',
    'order_tax', 'family', 'genus', 'species'
]

CULTURE_COLLECTION_COLUMNS = ['culture_collection']

HIGH_LEVEL_ENV_COLUMNS = [
    'animal_associated', 'aquatic', 'built', 'other',
    'sediment', 'soil', 'unassigned_high_level'
]

LOW_LEVEL_ENV_COLUMNS = [
    'human', 'invertebrate', 'other_vertebrate', 'unspecified_animal',
    'aquatic_other', 'aquatic_unspecified', 'freshwater', 'groundwater', 'marine',
    'built_other', 'drinking_water', 'wastewater',
    'air', 'bacteria', 'eukaryote_other', 'food', 'geothermal', 'hypersaline',
    'other_unspecified', 'plant_associated', 'subsurface', 'synthetic', 'viral',
    'freshwater_sediment', 'marine_sediment', 'sediment_unspecified',
    'desert', 'forest', 'rhizosphere', 'soil_agricultural', 'soil_other',
    'soil_unspecified', 'tundra_wetland', 'unassigned_low_level'
]

# ====================================
#
#   SQLite database setup
#
# ====================================
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

# ============================================================ #
# GFF parser                                                   #
# ============================================================ #
def extract_gene_info(line: list) -> tuple:
    """
    Extracts gene information from a GFF/GTF annotation line.

    Args:
        line (list): parsed annotation line with genomic feature data

    Returns:
        seqID (str): sequence identifier with underscores normalized
        annotation_ID (str): annotation identifier or 'NA' if not found
        parent (str): chromosome/contig name
        direction (str): strand direction (+/-)
        gene_start (str): start position
        gene_end (str): end position
        nuc_length (int): nucleotide sequence length
        aa_length (int): amino acid sequence length
    """
    # extract basic genomic coordinates and features
    parent = line[0]  # chromosome/contig name
    direction = line[6]  # strand direction (+/-)
    gene_start = line[3]  # start position
    gene_end = line[4]  # end position
    attributes = line[8]

    # parse attributes field
    attr_dict = {}
    for attr in attributes.split(';'):
        if '=' in attr:
            key, value = attr.split('=', 1)
            attr_dict[key] = value

    seqID = attr_dict.get('ID', '')
    annotation_ID = attr_dict.get('Name', '')


    # calculate sequence lengths
    nuc_length = abs(int(gene_end) - int(gene_start)) + 1  # nucleotide length
    aa_length = int(nuc_length / 3)  # amino acid length
    return seqID, annotation_ID, parent, direction, gene_start, gene_end, nuc_length, aa_length


def cugo_boundaries(direction: str,
                    prev_direction: str,
                    next_direction: str,
                    parent: str,
                    prev_parent: str,
                    next_parent: str,
                    cugo_count: int,
                    cugo_size: dict,
                    cugo_size_count: int,
                    prev_cugo: int = None):
    """
    Determines CUGO (Co-oriented Unidirectional Gene Order) boundaries and classifications.

    Args:
        direction (str): current gene strand direction (+/-)
        prev_direction (str): previous gene strand direction (+/-)
        next_direction (str): next gene strand direction (+/-)
        parent (str): current gene's chromosome/contig name
        prev_parent (str): previous gene's chromosome/contig name
        next_parent (str): next gene's chromosome/contig name
        cugo_count (int): running count of CUGO regions
        cugo_size (dict): mapping of CUGO IDs to their sizes
        cugo_size_count (int): running count of genes in current CUGO
        prev_cugo (int, optional): previous CUGO identifier

    Returns:
        cugo (int): current CUGO identifier
        cugo_start (str): start boundary type ('sequence_edge', 'strand_change', or 'NA')
        cugo_end (str): end boundary type ('sequence_edge', 'strand_change', or 'NA')
        cugo_count (int): updated CUGO count
        cugo_size_count (int): updated gene count for current CUGO
    """
    cugo, cugo_start, cugo_end = 'NA', 'NA', 'NA'

    # middle of cugo - same direction and parent for all three genes
    if (direction == prev_direction == next_direction) and (parent == prev_parent == next_parent):
        cugo = prev_cugo
        cugo_size_count += 1

    # start of contig/scaffold - different parent from previous gene
    elif parent != prev_parent:
        cugo = cugo_count
        cugo_size_count = 1

        if direction == '+':
            cugo_start = 'sequence_edge'
            if parent != next_parent:
                # single gene contig
                cugo_end = 'sequence_edge'
                cugo_size[cugo] = cugo_size_count
                cugo_size_count = 0
                cugo_count += 1
            elif direction != next_direction:
                # strand change at start of contig
                cugo_end = 'strand_change'
                cugo_size[cugo] = cugo_size_count
                cugo_size_count = 0
                cugo_count += 1
        else:  # direction == '-'
            cugo_end = 'sequence_edge'
            if parent != next_parent:
                # single gene contig
                cugo_start = 'sequence_edge'
                cugo_size[cugo] = cugo_size_count
                cugo_size_count = 0
                cugo_count += 1
            elif direction != next_direction:
                # strand change at start of contig
                cugo_start = 'strand_change'
                cugo_size[cugo] = cugo_size_count
                cugo_size_count = 0
                cugo_count += 1

    # end of contig/scaffold - different parent from next gene
    elif parent != next_parent:
        if direction == prev_direction:
            # continuing previous cugo
            cugo = prev_cugo
            cugo_size_count += 1
            if direction == '+':
                cugo_end = 'sequence_edge'
            else:
                cugo_start = 'sequence_edge'
            cugo_size[cugo] = cugo_size_count
            cugo_size_count = 0
            cugo_count += 1
        else:
            # new cugo at end of contig
            cugo = cugo_count
            cugo_size_count = 1
            if direction == '+':
                cugo_start = 'strand_change'
                cugo_end = 'sequence_edge'
            else:
                cugo_start = 'sequence_edge'
                cugo_end = 'strand_change'
            cugo_size[cugo] = cugo_size_count
            cugo_size_count = 0
            cugo_count += 1

    # strand change - different direction from previous or next gene
    elif direction != prev_direction or direction != next_direction:
        if direction != prev_direction and direction == next_direction:
            # start of new cugo due to strand change
            cugo = cugo_count
            cugo_size_count = 1
            if direction == '+':
                cugo_start = 'strand_change'
            else:
                cugo_end = 'strand_change'
        elif direction == prev_direction and direction != next_direction:
            # end of current cugo due to strand change
            cugo = prev_cugo
            cugo_size_count += 1
            if direction == '+':
                cugo_end = 'strand_change'
            else:
                cugo_start = 'strand_change'
            cugo_size[cugo] = cugo_size_count
            cugo_size_count = 0
            cugo_count += 1
        else:  # direction != prev_direction and direction != next_direction
            # isolated gene with different direction from both neighbors
            cugo = cugo_count
            cugo_size_count = 1
            cugo_start = cugo_end = 'strand_change'
            cugo_size[cugo] = cugo_size_count
            cugo_size_count = 0
            cugo_count += 1

    return cugo, cugo_start, cugo_end, cugo_count, cugo_size_count

def parse_single_gff(gff_content: str) -> list:
    """
    Parses GFF content and extracts gene information with CUGO boundary analysis.

    Args:
        gff_content (str): GFF format file content as string

    Returns:
        list: List of processed gene records, each containing:
            - seqID (str): gene/sequence identifier
            - parent_ID (str): chromosome/contig identifier
            - aa_length (int): amino acid sequence length
            - strand (str): gene strand direction (+/-)
            - annotation ID (str): functional category identifier
            - cugo_number (int): CUGO region identifier
    """
    reformat_data = []
    cugo_size = {}
    cugo_count = prev_direction = prev_parent = prev_feat_type = prev_cugo = 0
    cugo_size_count = 0
    prev_line = None
    cds_count = 0

    # split content into individual lines for processing
    lines = gff_content.strip().split('\n')

    for line_count, line in enumerate(lines, 1):
        # parse tab-delimited gff line
        clean_line = line.strip().split('\t')

        # skip malformed lines (gff requires 9 columns)
        if len(clean_line) != 9:
            continue

        # extract feature type and count cds features
        feat_type = clean_line[2]
        if feat_type == 'CDS':
            cds_count += 1

        # only process cds features or transitions from cds
        if feat_type != 'CDS' and prev_feat_type != 'CDS':
            continue
        prev_feat_type = feat_type

        # initialize with first cds line
        if prev_line is None:
            prev_line = clean_line
            continue

        # extract gene information from previous line
        seqID, annotation_ID, parent, direction, gene_start, gene_end, nuc_length, aa_length = extract_gene_info(
            prev_line)

        # get next gene's parent and direction for cugo boundary analysis
        next_parent = clean_line[0]
        next_direction = clean_line[6]

        # determine cugo boundaries and classifications
        cugo, cugo_start, cugo_end, cugo_count, cugo_size_count = cugo_boundaries(
            direction, prev_direction, next_direction,
            parent, prev_parent, next_parent,
            cugo_count, cugo_size, cugo_size_count, prev_cugo
        )

        # store processed gene data (all 11 fields for database)
        reformat_data.append([seqID, parent, aa_length,
                              direction, annotation_ID, cugo])

        # update tracking variables for next iteration
        prev_line = clean_line
        prev_direction = direction
        prev_parent = parent
        prev_cugo = cugo

    # process the last gene in the file (no next gene to compare)
    if len(clean_line) == 9 and clean_line[2] == 'CDS':
        seqID, annotation_ID, parent, direction, gene_start, gene_end, nuc_length, aa_length = extract_gene_info(
            clean_line)

        # handle last gene - no next gene, so use none for next_parent
        cugo, cugo_start, cugo_end, cugo_count, cugo_size_count = cugo_boundaries(
            direction, prev_direction, None,  # no next direction for last gene
            parent, prev_parent, None,  # no next parent for last gene
            cugo_count, cugo_size, cugo_size_count, prev_cugo
        )

        # store final gene data (all 6 fields for database)
        reformat_data.append([seqID, parent, aa_length,
                              direction, annotation_ID, cugo])

        # finalize cugo size tracking
        cugo_size[cugo] = cugo_size_count

    # Return complete data for database storage
    return reformat_data

def parse_gff_for_update(gff_content: str) -> list:
    """
    Parses GFF content for UPDATE operations (seqID and annotation_ID only).

    Args:
        gff_content (str): GFF format file content as string

    Returns:
        list: List of tuples (annotation_ID, seqID)
    """
    update_data = []
    lines = gff_content.strip().split('\n')

    for line in lines:
        clean_line = line.strip().split('\t')

        if len(clean_line) != 9 or clean_line[2] != 'CDS':
            continue

        attributes = clean_line[8]

        # parse attributes
        attr_dict = {}
        for attr in attributes.split(';'):
            if '=' in attr:
                key, value = attr.split('=', 1)
                attr_dict[key] = value

        seqID = attr_dict.get('ID', '')
        annotation_ID = attr_dict.get('Name', '')

        if seqID and annotation_ID:
            update_data.append((annotation_ID, seqID))

    return update_data



def process_gff_file(filepath: str,
                     update_mode: bool = False):
    """
    Process a single GFF file (gzipped or plain) from disk and parse it.

    Args:
        filepath (str): Full path to the GFF file.

    Returns:
        list: Parsed gene records.
    """
    try:
        filepath = Path(filepath)

        # Open gzipped or plain text
        if filepath.suffix == ".gz":
            with gzip.open(filepath, "rt") as f:
                content = f.read()
        else:
            with open(filepath, "r", encoding="utf-8") as f:
                content = f.read()

        if update_mode:
            return parse_gff_for_update(content)
        else:
            return parse_single_gff(content)
    except Exception as e:
        logger.warning(f"Failed to process {filepath}: {e}")
        return []


def process_tmhmm_file(tmhmm_filepath):
    try:
        with open(tmhmm_filepath, 'r') as f:
            tmhmm_data = []
            for line in f:
                if line.startswith('prot_ID'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    prot_id = parts[0]
                    try:
                        no_tmh = int(parts[1])
                    except ValueError:
                        no_tmh = 0
                    tmhmm_data.append((no_tmh, prot_id))
        return tmhmm_data
    except Exception:
        return []

def populate_taxonomy_table(conn: sqlite3.Connection,
                            taxonomy_filepath: str):
    try:
        filepath = Path(taxonomy_filepath)

        if filepath.suffix == '.gz':
            file_handle = gzip.open(filepath, 'rt')
        else:
            file_handle = open(filepath, 'r')

        with file_handle as f:
            next(f)

            taxonomy_data = []
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < len(TAXONOMY_COLUMNS):
                    continue

                row_values = tuple(parts[i] for i in range(len(TAXONOMY_COLUMNS)))

                taxonomy_data.append(row_values)

            placeholders = ', '.join(['?' for _ in TAXONOMY_COLUMNS])
            columns = ', '.join(TAXONOMY_COLUMNS)

            conn.executemany(f"""
                INSERT OR REPLACE INTO genome_data
                ({columns})
                VALUES ({placeholders}) 
            """, taxonomy_data)

        conn.commit()
        logger.info(f"Inserted taxonomy data for {len(taxonomy_data)} genomes")

    except Exception as e:
        logger.error(f"Error populating taxonomy data: {e}")

def populate_culture_collection(conn: sqlite3.Connection,
                                cc_file: str):
    try:
        with open(cc_file, 'r') as f:
            next(f)

            culture_collection_data = []
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 2:
                    continue

                genome_id = parts[0]
                culture_collection = int(parts[1])

                culture_collection_data.append((genome_id, culture_collection))

            query = """
                INSERT INTO genome_data (genome_ID, culture_collection)
                VALUES (?, ?)
                ON CONFLICT(genome_ID) DO UPDATE SET
                    culture_collection = excluded.culture_collection
            """

            conn.executemany(query, culture_collection_data)

            conn.commit()

    except Exception as e:
        logger.error(f"Error populating culture collection data: {e}")

def populate_high_level_environment(conn: sqlite3.Connection,
                                    high_level_environment_file: str):
    try:
        with open(high_level_environment_file, 'r') as f:
            next(f)

            high_level_environment_data = []
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < len(HIGH_LEVEL_ENV_COLUMNS) + 1:
                    continue

                genome_id = parts[0]

                env_values = [genome_id]
                for i, col in enumerate(HIGH_LEVEL_ENV_COLUMNS, start=1):
                    value = float(parts[i]) if parts[i].strip() else None
                    env_values.append(value)

                high_level_environment_data.append(tuple(env_values))

            columns = ['genome_ID'] + HIGH_LEVEL_ENV_COLUMNS
            placeholders = ', '.join(['?' for _ in columns])
            column_names = ', '.join(columns)
            update_clause = ', '.join([f'{col} = excluded.{col}' for col in HIGH_LEVEL_ENV_COLUMNS])

            query = f"""
                INSERT INTO genome_data ({column_names})
                VALUES ({placeholders})
                ON CONFLICT(genome_ID) DO UPDATE SET
                    {update_clause}
            """

            conn.executemany(query, high_level_environment_data)
            conn.commit()
            logger.info(f"Processed high level environment data for {len(high_level_environment_data)} genomes")

    except Exception as e:
        logger.error(f"Error populating high level environment data: {e}")


def populate_low_level_environment(conn: sqlite3.Connection,
                                   low_level_environment_file: str):
    try:
        with open(low_level_environment_file, 'r') as f:
            next(f)

            low_level_environment_data = []
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < len(LOW_LEVEL_ENV_COLUMNS) + 1:
                    continue

                genome_id = parts[0]

                env_values = [genome_id]
                for i, col in enumerate(LOW_LEVEL_ENV_COLUMNS, start=1):
                    value = float(parts[i]) if parts[i].strip() else None
                    env_values.append(value)

                low_level_environment_data.append(tuple(env_values))

        columns = ['genome_ID'] + LOW_LEVEL_ENV_COLUMNS
        placeholders = ', '.join(['?' for _ in columns])
        column_names = ', '.join(columns)
        update_clause = ', '.join([f'{col} = excluded.{col}' for col in LOW_LEVEL_ENV_COLUMNS])

        query = f"""
            INSERT INTO genome_data ({column_names})
            VALUES ({placeholders})
            ON CONFLICT(genome_ID) DO UPDATE SET
                {update_clause}
        """

        conn.executemany(query, low_level_environment_data)
        conn.commit()
        logger.info(f"Processed low level environment data for {len(low_level_environment_data)} genomes")

    except Exception as e:
        logger.error(f"Error populating low level environment data: {e}")

def stream_all_proteins(fasta: str,
                        batch_size: int = 10000):
    try:
        filepath = Path(fasta)

        if filepath.suffix == '.gz':
            file_handle = gzip.open(filepath, 'rt')
        else:
            file_handle = open(filepath, 'r')

        batch = []
        current_id = None

        with file_handle as f:
            for line in f:
                line = line.strip()

                if not line:
                    continue

                if line.startswith('>'):
                    current_id = line[1:].split()[0]
                else:
                    if current_id is not None:
                        compressed = compress_sequence(line)
                        batch.append((current_id, compressed))
                        current_id = None

                        if len(batch) >= batch_size:
                            yield batch
                            batch = []

            if batch:
                yield batch

    except Exception as e:
        logger.warning(f"Failed to process protein FASTA {fasta}: {e}")

def populate_protein_sequences(conn: sqlite3.Connection,
                               protein_fasta_path: str,
                               ):
    logger.info(f"Processing protein FASTA file: {protein_fasta_path}")

    total_sequences = 0

    for batch in stream_all_proteins(protein_fasta_path):
        insert_data = [(seqid, seq_blob) for seqid, seq_blob in batch]
        conn.executemany("""
                         INSERT INTO protein_data (seqID, protein_seq)
                         VALUES (?, ?)
                             ON CONFLICT(seqID) DO UPDATE SET
                             protein_seq = excluded.protein_seq
                         """, insert_data)
        conn.commit()
        total_sequences += len(batch)

    logger.info(f"Processed {total_sequences} protein sequences")

def database(cog_gff_tar_path: str,
          kegg_gff_tar_path: str = None,
          pfam_gff_tar_path: str = None,
          tmhmm_tar_path: str = None,
          protein_fasta_path: str = None,
          taxonomy_path: str = None,
          culture_collection_path: str = None,
          high_level_environment_path: str = None,
          low_level_environment_path: str = None,
          output_dir: str = None,
          globdb_version: str = None,
          ):
    """
    Main parsing function - processes COG first with full data, then updates with KEGG and Pfam.
    All annotations use Name= attribute from GFF files.
    """
    from logging import getLogger
    logger = getLogger(__name__)

    logger.info("Processing files: COG → TMHMM → KEGG → Pfam → Protein sequences → Taxonomy → Culture collection "
                "→ High level environment data → Low level environment data ")

    # db_path = ensure_path(target=f"globdb_{globdb_version}_cugo.db")
    db_path = f"globdb_{globdb_version}_cugo.db"

    if output_dir is None:
        output_dir = '.'

    output_path = Path(output_dir)
    tempdir = output_path / 'tempdir'
    if tempdir.exists():
        shutil.rmtree(tempdir)
    tempdir.mkdir(parents=True, exist_ok=True)

    # Setup database
    conn = setup_database(db_path)

    # ===== STEP 1: Process COG GFF files =====
    logger.info("STEP 1: Processing COG GFF files...")
    subprocess.run(["tar", "-xzf", cog_gff_tar_path, "-C", str(tempdir)], check=True)

    cog_gff_files = list(tempdir.rglob("*.gff.gz")) + list(tempdir.rglob("*.gff"))
    logger.info(f"Found {len(cog_gff_files)} COG GFF files")

    for gff_file in tqdm(cog_gff_files, desc="Processing COG GFF"):
        gff_data = process_gff_file(str(gff_file), update_mode=False)
        if gff_data:
            conn.executemany("""
                INSERT OR REPLACE INTO protein_data
                (seqID, parent_ID, aa_length, strand, COG_ID, cugo_number)
                VALUES (?, ?, ?, ?, ?, ?)
            """, gff_data)
            conn.commit()

    logger.info("Creating protein_data indexes...")
    conn.execute('CREATE INDEX IF NOT EXISTS idx_all_seqid ON protein_data(seqID)')
    conn.execute('CREATE INDEX IF NOT EXISTS idx_parent_id ON protein_data(parent_ID)')
    conn.commit()
    logger.info("protein_data indexes created")

    shutil.rmtree(tempdir)
    tempdir.mkdir(parents=True, exist_ok=True)

    # ===== STEP 2: Process TMHMM files =====
    if tmhmm_tar_path:
        logger.info("STEP 2: Processing TMHMM files...")
        subprocess.run(["tar", "-xzf", tmhmm_tar_path, "-C", str(tempdir)], check=True)

        tmhmm_files = list(tempdir.rglob("*_tmhmm_clean"))
        logger.info(f"Found {len(tmhmm_files)} TMHMM files")

        for tmhmm_file in tqdm(tmhmm_files, desc="Processing TMHMM"):
            tmhmm_data = process_tmhmm_file(str(tmhmm_file))
            if tmhmm_data:
                # Changed to handle new entries
                conn.executemany("""
                                 INSERT INTO protein_data (seqID, no_tmh)
                                 VALUES (?, ?)
                                     ON CONFLICT(seqID) DO UPDATE SET
                                     no_tmh = excluded.no_tmh
                                 """, [(seqid, no_tmh) for no_tmh, seqid in tmhmm_data])
                conn.commit()

        shutil.rmtree(tempdir)
        tempdir.mkdir(parents=True, exist_ok=True)

    # ===== STEP 3: Process KEGG GFF files =====
    if kegg_gff_tar_path:
        logger.info("STEP 3: Updating KEGG annotations...")
        subprocess.run(["tar", "-xzf", kegg_gff_tar_path, "-C", str(tempdir)], check=True)

        kegg_gff_files = list(tempdir.rglob("*.gff.gz")) + list(tempdir.rglob("*.gff"))
        logger.info(f"Found {len(kegg_gff_files)} KEGG GFF files")

        for gff_file in tqdm(kegg_gff_files, desc="Updating KEGG IDs"):
            kegg_data = process_gff_file(str(gff_file), update_mode=True)
            if kegg_data:
                conn.executemany("""
                                 INSERT INTO protein_data (seqID, KEGG_ID)
                                 VALUES (?, ?)
                                     ON CONFLICT(seqID) DO UPDATE SET
                                     KEGG_ID = excluded.KEGG_ID
                                 """, [(seqid, kegg_id) for kegg_id, seqid in kegg_data])
                conn.commit()

        shutil.rmtree(tempdir)
        tempdir.mkdir(parents=True, exist_ok=True)

    # ===== STEP 4: Process Pfam GFF files =====
    if pfam_gff_tar_path:
        logger.info("STEP 4: Updating Pfam annotations...")
        subprocess.run(["tar", "-xzf", pfam_gff_tar_path, "-C", str(tempdir)], check=True)

        pfam_gff_files = list(tempdir.rglob("*.gff.gz")) + list(tempdir.rglob("*.gff"))
        logger.info(f"Found {len(pfam_gff_files)} Pfam GFF files")

        for gff_file in tqdm(pfam_gff_files, desc="Updating Pfam IDs"):
            pfam_data = process_gff_file(str(gff_file), update_mode=True)
            if pfam_data:
                conn.executemany("""
                                 INSERT INTO protein_data (seqID, Pfam_ID)
                                 VALUES (?, ?)
                                     ON CONFLICT(seqID) DO UPDATE SET
                                     Pfam_ID = excluded.Pfam_ID
                                 """, [(seqid, pfam_id) for pfam_id, seqid in pfam_data])
                conn.commit()

        shutil.rmtree(tempdir)

    # ===== STEP 5: Process Protein FASTA file =====
    if protein_fasta_path:
        logger.info("STEP 5: Processing protein sequences...")
        populate_protein_sequences(conn, protein_fasta_path)

    # ===== STEP 6: Process Taxonomy file =====
    if taxonomy_path:
        logger.info("STEP 6: Processing taxonomy data...")
        populate_taxonomy_table(conn, taxonomy_path)

    logger.info("Creating genome_data index...")
    conn.execute('CREATE INDEX IF NOT EXISTS idx_genome_id ON genome_data(genome_ID)')
    conn.commit()
    logger.info("genome_data index created")

    # ==== STEP 7: Process culture collection file ====
    if culture_collection_path:
        logger.info("STEP 7: Processing culture collection data...")
        populate_culture_collection(conn, culture_collection_path)

    # ==== STEP 8: Process high level environment data ====
    if high_level_environment_path:
        logger.info("STEP 8: Processing high level environment data...")
        populate_high_level_environment(conn, high_level_environment_path)

    # ==== STEP 9: Process low level environment data ====
    if low_level_environment_path:
        logger.info("STEP 9: Processing low level environment data...")
        populate_low_level_environment(conn, low_level_environment_path)

    logger.info("Optimizing database with VACUUM...")
    conn.execute("VACUUM")
    conn.commit()
    logger.info("VACUUM complete")

    conn.close()
    logger.info("Database creation complete!")


# ======================================
# METADATA RETRIEVAL
# ======================================
def build_query(include_taxonomy: bool = False,
                include_annotation: bool = False,
                include_culture_collection: bool = False,
                include_high_level_environment: bool = False,
                include_low_level_environment: bool = False,
                batch_size: int = 1) -> str:
    select_cols = [f'p.{col}' for col in BASE_COLUMNS]

    if include_annotation:
        select_cols.extend([f'p.{col}' for col in ANNOTATION_COLUMNS])

    if include_taxonomy:
        select_cols.extend([f'g.{col}' for col in TAXONOMY_COLUMNS])

    if include_culture_collection:
        select_cols.extend([f'g.{col}' for col in CULTURE_COLLECTION_COLUMNS])

    if include_high_level_environment:
        select_cols.extend([f'g.{col}' for col in HIGH_LEVEL_ENV_COLUMNS])

    if include_low_level_environment:
        select_cols.extend([f'g.{col}' for col in LOW_LEVEL_ENV_COLUMNS])

    query = f"SELECT {', '.join(select_cols)} FROM protein_data p"

    needs_genome_join = (
            include_taxonomy or include_culture_collection or
            include_high_level_environment or include_low_level_environment
    )

    if needs_genome_join:
        query += """
            LEFT JOIN genome_data g ON 
            CASE 
                WHEN instr(p.seqID, '___') > 0 
                THEN substr(p.seqID, 1, instr(p.seqID, '___') - 1)
                ELSE p.seqID
            END = g.genome_ID
            """

    if batch_size > 1:
        placeholders = ','.join(['?' for _ in range(batch_size)])
        query += f" WHERE p.seqID IN ({placeholders})"
    else:
        query += " WHERE p.seqID = ?"

    return query


def get_header(include_taxonomy: bool,
               include_annotation: bool,
               include_culture_collection: bool,
               include_high_level_environment: bool,
               include_low_level_environment: bool):

    header = BASE_COLUMNS.copy()

    if include_annotation:
        header.extend(ANNOTATION_COLUMNS)

    if include_taxonomy:
        header.extend(TAXONOMY_COLUMNS)

    if include_culture_collection:
        header.extend(CULTURE_COLLECTION_COLUMNS)

    if include_high_level_environment:
        header.extend(HIGH_LEVEL_ENV_COLUMNS)

    if include_low_level_environment:
        header.extend(LOW_LEVEL_ENV_COLUMNS)

    return header


def process_batch(db_path: str,
                  batch: list,
                  include_taxonomy: bool,
                  include_annotation: bool,
                  include_culture_collection: bool,
                  include_high_level_environment: bool,
                  include_low_level_environment: bool):

    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    query = build_query(
        include_taxonomy=include_taxonomy,
        include_annotation=include_annotation,
        include_culture_collection=include_culture_collection,
        include_high_level_environment=include_high_level_environment,
        include_low_level_environment=include_low_level_environment,
        batch_size=len(batch)
    )

    cursor.execute(query, batch)
    results = cursor.fetchall()

    result_dict = {row[0]: row for row in results}

    conn.close()
    return batch, result_dict


def meta(db_path: str,
         fasta: str,
         output: str,
         threads: int = 1,
         include_annotation: bool = False,
         include_taxonomy: bool = False,
         include_culture_collection: bool = False,
         include_high_level_environment: bool = False,
         include_low_level_environment: bool = False,
         all_metadata: bool = False,
         force: bool = False):

    if fasta:
        seq_dict = read_fasta_to_dict(fasta)
        seq_ids = list(seq_dict.keys())
        dataset_name = determine_dataset_name(fasta, '.', 0)
        output_path = ensure_path(
            output,
            f'{dataset_name}_metadata.tsv',
            force=force
        )
    else:
        logger.error(
            "No valid input file found. Please specify path to protein FASTA"
        )
        return

    # Enable everything if --all-metadata is set
    if all_metadata:
        include_annotation = True
        include_taxonomy = True
        include_culture_collection = True
        include_high_level_environment = True
        include_low_level_environment = True

    # Ensure at least one metadata type is requested
    if not any([
        include_annotation,
        include_taxonomy,
        include_culture_collection,
        include_high_level_environment,
        include_low_level_environment
    ]):
        logger.error(
            "Must specify at least one metadata type "
            "(--taxonomy, --annotation, --culture-collection, "
            "--high-level-environment, --low-level-environment, or --all-metadata)"
        )
        return

    BATCH_SIZE = 500

    header = get_header(
        include_taxonomy,
        include_annotation,
        include_culture_collection,
        include_high_level_environment,
        include_low_level_environment
    )

    # Split seq_ids into batches
    batches = [
        seq_ids[i:i + BATCH_SIZE]
        for i in range(0, len(seq_ids), BATCH_SIZE)
    ]

    found_ids = set()
    missing_ids = set()
    batch_results = {}

    # Process batches in parallel
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {
            executor.submit(
                process_batch,
                db_path,
                batch,
                include_taxonomy,
                include_annotation,
                include_culture_collection,
                include_high_level_environment,
                include_low_level_environment
            ): idx
            for idx, batch in enumerate(batches)
        }

        for future in tqdm(
                as_completed(futures),
                total=len(batches),
                desc="Querying database"
        ):
            batch_idx = futures[future]
            try:
                batch, result_dict = future.result()
                batch_results[batch_idx] = (batch, result_dict)
            except Exception as e:
                logger.error(
                    f"Error processing batch {batch_idx}: {e}"
                )

    # Write results in original order
    with open(output_path, 'w') as out:
        out.write('\t'.join(header) + '\n')

        for batch_idx in sorted(batch_results.keys()):
            batch, result_dict = batch_results[batch_idx]

            for seq_id in batch:
                if seq_id in result_dict:
                    found_ids.add(seq_id)
                    row = [
                        str(x) if x is not None else ''
                        for x in result_dict[seq_id]
                    ]
                    out.write('\t'.join(row) + '\n')
                else:
                    missing_ids.add(seq_id)

    logger.info(f"Metadata has been saved to {output_path}")

    if missing_ids:
        logger.warning(
            f"{len(missing_ids)} sequence IDs were not found in the database"
        )



def metadata_categories():
    print("───────────────────────────")
    print("Protein metadata")
    print("───────────────────────────")
    print("  • sequence information")
    for item in BASE_COLUMNS[1:]:
        print(f"      - {item}")
    print()
    print("  • annotations")
    for item in ANNOTATION_COLUMNS:
        print(f"      - {item}")
    print()  # blank line

    print("───────────────────────────")
    print("Genome metadata")
    print("───────────────────────────")
    print("  • taxonomy")
    for item in TAXONOMY_COLUMNS:
        print(f"      - {item}")
    print()

    print("  • culture information")
    for item in CULTURE_COLLECTION_COLUMNS:
        print(f"      - {item}")
    print()

    print("  • high-level environment data")
    for item in HIGH_LEVEL_ENV_COLUMNS:
        print(f"      - {item}")
    print()

    print("  • low-level environment data")
    for item in LOW_LEVEL_ENV_COLUMNS:
        print(f"      - {item}")
    print()