from .cugo import *
from .util import *

import sqlite3
import subprocess
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed

logger = logging.getLogger(__name__)

# ====================================
#
#   SQLite database setup
#
# ====================================
def setup_database(db_path: str) -> sqlite3.Connection:
    conn = sqlite3.connect(db_path)

    conn.execute('''
                 CREATE TABLE IF NOT EXISTS protein_data (
                     seqID TEXT PRIMARY KEY,
                     parent_ID TEXT,
                     aa_length INTEGER,
                     strand TEXT,
                     COG_ID TEXT,
                     KEGG_ID TEXT,
                     Pfam_ID TEXT,
                     cugo_number INTEGER,
                     no_tmh INTEGER,
                     protein_seq BLOB
                 )
                 ''')

    conn.execute('''
                 CREATE TABLE IF NOT EXISTS genome_data (
                     genome_ID TEXT PRIMARY KEY,
                     domain TEXT,
                     phylum TEXT,
                     class TEXT,
                     order_tax TEXT,
                     family TEXT,
                     genus TEXT,
                     species TEXT,
                     culture_collection INT,
                     animal_associated REAL,
                     aquatic REAL,
                     built REAL,
                     other REAL,
                     sediment REAL,
                     soil REAL,
                     unassigned_high_level REAL,
                     human REAL,
                     invertebrate REAL,
                     other_vertebrate REAL,
                     unspecified_animal REAL,
                     aquatic_other REAL,
                     aquatic_unspecified REAL,
                     freshwater REAL,
                     groundwater REAL,
                     marine REAL,
                     built_other REAL,
                     drinking_water REAL,
                     wastewater REAL,
                     air REAL,
                     bacteria REAL,
                     eukaryote_other REAL,
                     food REAL,
                     geothermal REAL,
                     hypersaline REAL,
                     other_unspecified REAL,
                     plant_associated REAL,
                     subsurface REAL,
                     synthetic REAL,
                     viral REAL,
                     freshwater_sediment REAL,
                     marine_sediment REAL,
                     sediment_unspecified REAL,
                     desert REAL,
                     forest REAL,
                     rhizosphere REAL,
                     soil_agricultural REAL,
                     soil_other REAL,
                     soil_unspecified REAL,
                     tundra_wetland REAL,
                     unassigned_low_level REAL     
                 )
                 ''')

    conn.execute('''
                 CREATE VIEW IF NOT EXISTS protein_data_readable AS
                 SELECT
                     seqID,
                     parent_ID,
                     aa_length,
                     strand,
                     COG_ID,
                     KEGG_ID,
                     Pfam_ID,
                     cugo_number,
                     no_tmh,
                     CASE WHEN protein_seq IS NULL THEN NULL ELSE '<COMPRESSED_BLOB>' END as protein_seq_status
                 FROM protein_data
                 ''')


    conn.commit()
    return conn

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
                if len(parts) < 8:
                    continue

                genome_id = parts[0]
                domain = parts[1]
                phylum = parts[2]
                class_ = parts[3]
                order = parts[4]
                family = parts[5]
                genus = parts[6]
                species = parts[7]

                taxonomy_data.append((genome_id, domain, phylum, class_, order, family, genus, species))

            conn.executemany("""
                INSERT OR REPLACE INTO genome_data
                (genome_id, domain, phylum, class, order_tax, family, genus, species)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?) 
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

                culture_collection_data.append((culture_collection, genome_id))

            conn.executemany("""
                UPDATE genome_data
                SET culture_collection = ?
                WHERE genome_id = ?
            """, culture_collection_data)
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
                if len(parts) < 9:
                    continue

                genome_id = parts[0]
                animal_associated = float(parts[1]) if parts[1].strip() else None
                aquatic = float(parts[2]) if parts[2].strip() else None
                built = float(parts[3]) if parts[3].strip() else None
                other = float(parts[4]) if parts[4].strip() else None
                sediment = float(parts[5]) if parts[5].strip() else None
                soil = float(parts[6]) if parts[6].strip() else None
                unassigned_high_level = float(parts[7]) if parts[7].strip() else None

                high_level_environment_data.append((unassigned_high_level, soil, sediment, other, built, aquatic, animal_associated, genome_id))

            conn.executemany("""
                UPDATE genome_data
                SET 
                    unassigned_high_level = ?,
                    soil = ?,
                    sediment = ?,
                    other = ?,
                    built = ?,
                    aquatic = ?,
                    animal_associated = ?
                WHERE genome_id = ?
            """, high_level_environment_data)
            conn.commit()

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
                if len(parts) < 36:
                    continue

                genome_id = parts[0]
                human = float(parts[1]) if parts[1].strip() else None
                invertebrate = float(parts[2]) if parts[2].strip() else None
                other_vertebrate = float(parts[3]) if parts[3].strip() else None
                unspecified_animal = float(parts[4]) if parts[4].strip() else None
                aquatic_other = float(parts[5]) if parts[5].strip() else None
                aquatic_unspecified = float(parts[6]) if parts[6].strip() else None
                freshwater = float(parts[7]) if parts[7].strip() else None
                groundwater = float(parts[8]) if parts[8].strip() else None
                marine = float(parts[9]) if parts[9].strip() else None
                built_other = float(parts[10]) if parts[10].strip() else None
                drinking_water = float(parts[11]) if parts[11].strip() else None
                wastewater = float(parts[12]) if parts[12].strip() else None
                air = float(parts[13]) if parts[13].strip() else None
                bacteria = float(parts[14]) if parts[14].strip() else None
                eukaryote_other = float(parts[15]) if parts[15].strip() else None
                food = float(parts[16]) if parts[16].strip() else None
                geothermal = float(parts[17]) if parts[17].strip() else None
                hypersaline = float(parts[18]) if parts[18].strip() else None
                other_unspecified = float(parts[19]) if parts[19].strip() else None
                plant_associated = float(parts[20]) if parts[20].strip() else None
                subsurface = float(parts[21]) if parts[21].strip() else None
                synthetic = float(parts[22]) if parts[22].strip() else None
                viral = float(parts[23]) if parts[23].strip() else None
                freshwater_sediment = float(parts[24]) if parts[24].strip() else None
                marine_sediment = float(parts[25]) if parts[25].strip() else None
                sediment_unspecified = float(parts[26]) if parts[26].strip() else None
                desert = float(parts[27]) if parts[27].strip() else None
                forest = float(parts[28]) if parts[28].strip() else None
                rhizosphere = float(parts[29]) if parts[29].strip() else None
                soil_agricultural = float(parts[30]) if parts[30].strip() else None
                soil_other = float(parts[31]) if parts[31].strip() else None
                soil_unspecified = float(parts[32]) if parts[32].strip() else None
                tundra_wetland = float(parts[33]) if parts[33].strip() else None
                unassigned_low_level = float(parts[34]) if parts[34].strip() else None

                low_level_environment_data.append((unassigned_low_level, tundra_wetland, soil_unspecified, soil_other,
                                                   soil_agricultural, rhizosphere, forest, desert, sediment_unspecified,
                                                   marine_sediment, freshwater_sediment, viral, synthetic, subsurface,
                                                   plant_associated, other_unspecified, hypersaline, geothermal, food,
                                                   eukaryote_other,bacteria, air, wastewater, drinking_water, built_other,
                                                   marine, groundwater, freshwater, aquatic_unspecified, aquatic_other,
                                                   unspecified_animal, other_vertebrate, invertebrate, human, genome_id
                                                   ))

            conn.executemany("""
                UPDATE genome_data
                SET
                    unassigned_low_level = ?,
                    tundra_wetland = ?, 
                    soil_unspecified = ?, 
                    soil_other = ?,
                    soil_agricultural = ?,
                    rhizosphere = ?,
                    forest = ?,
                    desert = ?, 
                    sediment_unspecified = ?,
                    marine_sediment = ?, 
                    freshwater_sediment = ?, 
                    viral = ?, 
                    synthetic = ?, 
                    subsurface = ?,
                    plant_associated = ?, 
                    other_unspecified = ?, 
                    hypersaline = ?, 
                    geothermal = ?, 
                    food = ?,
                    eukaryote_other = ?, 
                    bacteria = ?, 
                    air = ?, 
                    wastewater = ?, 
                    drinking_water = ?, 
                    built_other = ?,
                    marine = ?,
                    groundwater = ?, 
                    freshwater = ?, 
                    aquatic_unspecified = ?, 
                    aquatic_other = ?,
                    unspecified_animal = ?, 
                    other_vertebrate = ?, 
                    invertebrate = ?, 
                    human = ?
                WHERE genome_id = ?
            """, low_level_environment_data)
            conn.commit()

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
        update_data = [(seq_blob, seqid) for seqid, seq_blob in batch]
        conn.executemany("""
                         UPDATE protein_data
                         SET protein_seq = ?
                         WHERE seqID = ?
                         """, update_data)
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
                conn.executemany("""
                                 UPDATE protein_data
                                 SET no_tmh = ?
                                 WHERE seqID = ?
                                 """, tmhmm_data)
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
                                 UPDATE protein_data
                                 SET KEGG_ID = ?
                                 WHERE seqID = ?
                                 """, kegg_data)
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
                                 UPDATE protein_data
                                 SET Pfam_ID = ?
                                 WHERE seqID = ?
                                 """, pfam_data)
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
def build_query(include_taxonomy: bool, include_annotation: bool, batch_size: int = 1) -> str:
    select_cols = ['p.seqID', 'p.parent_ID', 'p.aa_length', 'p.strand']

    if include_annotation:
        select_cols.extend(['p.COG_ID', 'p.KEGG_ID', 'p.Pfam_ID'])

    if include_taxonomy:
        select_cols.extend([
            'g.genome_ID', 'g.domain', 'g.phylum', 'g.class',
            'g.order_tax', 'g.family', 'g.genus', 'g.species'
        ])

    query = f"SELECT {', '.join(select_cols)} FROM protein_data p"

    if include_taxonomy:
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

def get_header(include_taxonomy: bool, include_annotation: bool):
    header = ['seqID', 'parent_ID', 'aa_length', 'strand']

    if include_annotation:
        header.extend(['COG_ID', 'KEGG_ID', 'Pfam_ID'])

    if include_taxonomy:
        header.extend(['genome_ID', 'domain', 'phylum', 'class',
                       'order', 'family', 'genus', 'species'])

    return header

def process_batch(db_path: str,
                  batch: list,
                  include_taxonomy: bool,
                  include_annotation: bool):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    query = build_query(include_taxonomy, include_annotation, len(batch))
    cursor.execute(query, batch)
    results = cursor.fetchall()

    result_dict = {}
    for result in results:
        seq_id = result[0]
        result_dict[seq_id] = result

    conn.close()
    return batch, result_dict

def meta(db_path: str,
             fasta: str,
             output: str,
             threads: int = 1,
             include_annotation: bool = False,
             include_taxonomy: bool = False,
             all_metadata: bool = False,
             force: bool = False):
    if fasta:
        seq_dict = read_fasta_to_dict(fasta)
        seq_ids = list(seq_dict.keys())
        dataset_name = determine_dataset_name(fasta, '.', 0)
        output_path = ensure_path(output, f'{dataset_name}_metadata.tsv', force=force)
    else:
        logger.error("No valid input file found. Please specify path to protein FASTA")
        return

    if all_metadata:
        include_annotation = True
        include_taxonomy = True

    if not (include_annotation or include_taxonomy):
        logger.error("Must specify at least one data type (--taxonomy, --annotation, or --all)")
        return

    conn = sqlite3.connect(db_path)

    BATCH_SIZE = 500

    header = get_header(include_taxonomy, include_annotation)

    # Split seq_ids into batches
    batches = [seq_ids[i:i + BATCH_SIZE] for i in range(0, len(seq_ids), BATCH_SIZE)]

    found_ids = set()
    missing_ids = set()

    # Store results with batch index to maintain order
    batch_results = {}

    # Process batches in parallel
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(process_batch, db_path, batch, include_taxonomy, include_annotation): idx
                   for idx, batch in enumerate(batches)}

        for future in tqdm(as_completed(futures), total=len(batches), desc="Querying database"):
            batch_idx = futures[future]
            try:
                batch, result_dict = future.result()
                batch_results[batch_idx] = (batch, result_dict)
            except Exception as e:
                logger.error(f"Error processing batch {batch_idx}: {e}")

    # Write results in original order
    with open(output_path, 'w') as out:
        out.write('\t'.join(header) + '\n')

        for batch_idx in sorted(batch_results.keys()):
            batch, result_dict = batch_results[batch_idx]

            for seq_id in batch:
                if seq_id in result_dict:
                    found_ids.add(seq_id)
                    row = [str(x) if x is not None else '' for x in result_dict[seq_id]]
                    out.write('\t'.join(row) + '\n')
                else:
                    missing_ids.add(seq_id)

        logger.info(f"Metadata has been saved to {output_path}")

    conn.close()