#!/usr/bin/env python3

from .util import *

import pandas as pd
import logging
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable, get_cmap
from typing import Optional
from pathlib import Path
import yaml
import tarfile
import gzip
import sqlite3
import multiprocessing as mp
from tqdm import tqdm
import subprocess

logger = logging.getLogger(__name__)

# ===========================
# CUGO GFF PARSER
# ===========================
def setup_database(db_path: str) -> sqlite3.Connection:
    conn = sqlite3.connect(db_path)

    conn.execute('''
        CREATE TABLE IF NOT EXISTS gff_data (
            seqID TEXT PRIMARY KEY,
            parent_ID TEXT,
            gene_start INTEGER,
            gene_end INTEGER,
            nuc_length INTEGER,
            aa_length INTEGER,
            strand TEXT,
            COG_ID TEXT,
            cugo_number INTEGER,
            cugo_start TEXT,
            cugo_end TEXT
        )
    ''')

    conn.execute('''
            CREATE TABLE IF NOT EXISTS tmhmm_data (
                seqID TEXT PRIMARY KEY,
                no_tmh INTEGER
            )
        ''')

    conn.execute('CREATE INDEX IF NOT EXISTS idx_gff_seqid ON gff_data(seqID)')
    conn.execute('CREATE INDEX IF NOT EXISTS idx_tmhmm_seqid ON tmhmm_data(seqID)')

    conn.commit()
    return conn

def extract_gene_info(line: list) -> tuple:
    """
    Extracts gene information from a GFF/GTF annotation line.

    Args:
        line (list): parsed annotation line with genomic feature data
        tmhmm_dict (dict): mapping of protein IDs to transmembrane helix counts

    Returns:
        seqID (str): sequence identifier with underscores normalized
        COG_ID (str): COG identifier or 'NA' if not found
        parent (str): chromosome/contig name
        direction (str): strand direction (+/-)
        gene_start (str): start position
        gene_end (str): end position
        nuc_length (int): nucleotide sequence length
        aa_length (int): amino acid sequence length
    """
    # parse annotation field (9th column) for sequence and COG identifiers
    annotation = line[8].split(';')
    seqID = annotation[0].split('=')[1].replace('__', '_')
    COG_ID = annotation[1].split('=')[1] if len(annotation) > 1 else 'NA'

    # extract basic genomic coordinates and features
    parent = line[0]  # chromosome/contig name
    direction = line[6]  # strand direction (+/-)
    gene_start = line[3]  # start position
    gene_end = line[4]  # end position

    # calculate sequence lengths
    nuc_length = abs(int(gene_end) - int(gene_start)) + 1  # nucleotide length
    aa_length = int(nuc_length / 3)  # amino acid length
    return seqID, COG_ID, parent, direction, gene_start, gene_end, nuc_length, aa_length


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
            - COG_ID (str): COG functional category identifier
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
        seqID, COG_ID, parent, direction, gene_start, gene_end, nuc_length, aa_length = extract_gene_info(
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
        reformat_data.append([seqID, parent, gene_start, gene_end, nuc_length, aa_length,
                              direction, COG_ID, cugo, cugo_start, cugo_end])

        # update tracking variables for next iteration
        prev_line = clean_line
        prev_direction = direction
        prev_parent = parent
        prev_cugo = cugo

    # process the last gene in the file (no next gene to compare)
    if len(clean_line) == 9 and clean_line[2] == 'CDS':
        seqID, COG_ID, parent, direction, gene_start, gene_end, nuc_length, aa_length = extract_gene_info(
            clean_line)

        # handle last gene - no next gene, so use none for next_parent
        cugo, cugo_start, cugo_end, cugo_count, cugo_size_count = cugo_boundaries(
            direction, prev_direction, None,  # no next direction for last gene
            parent, prev_parent, None,  # no next parent for last gene
            cugo_count, cugo_size, cugo_size_count, prev_cugo
        )

        # store final gene data (all 11 fields for database)
        reformat_data.append([seqID, parent, gene_start, gene_end, nuc_length, aa_length,
                              direction, COG_ID, cugo, cugo_start, cugo_end])

        # finalize cugo size tracking
        cugo_size[cugo] = cugo_size_count

    # Return complete data for database storage
    return reformat_data


def process_gff_file(filepath: str):
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

        gff_data = parse_single_gff(content)
        return gff_data
    except Exception as e:
        # Optional: log the error
        # logger.warning(f"Failed to process {filepath}: {e}")
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
                    prot_id = parts[0].replace('___', '__')
                    try:
                        no_tmh = int(parts[1])
                    except ValueError:
                        no_tmh = 0
                    tmhmm_data.append((prot_id, no_tmh))
        return tmhmm_data
    except Exception:
        return []

def export_to_tsv(db_path: str, output_tsv: str):
    """Export merged data from SQLite to TSV"""
    conn = sqlite3.connect(db_path)

    # Join GFF and TMHMM data
    query = '''
        SELECT 
            g.seqID,
            g.parent_ID,
            g.aa_length,
            g.strand,
            g.COG_ID,
            g.cugo_number,
            COALESCE(t.no_tmh, 0) as no_TMH
        FROM gff_data g
        LEFT JOIN tmhmm_data t ON g.seqID = t.seqID
        ORDER BY g.seqID
    '''

    cursor = conn.execute(query)

    # Write TSV
    with open(output_tsv, 'w') as f:
        # Header
        f.write('seqID\tparent_ID\taa_length\tstrand\tCOG_ID\tCUGO_number\tno_TMH\n')

        # Data
        for row in cursor:
            f.write('\t'.join(map(str, row)) + '\n')

    conn.close()


# !/usr/bin/env python3

from .util import *

import pandas as pd
import logging
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable, get_cmap
from typing import Optional
from pathlib import Path
import yaml
import tarfile
import gzip
import sqlite3
from tqdm import tqdm
import subprocess
import shutil

logger = logging.getLogger(__name__)


# ===========================
# CUGO GFF PARSER
# ===========================
def setup_database(db_path: str) -> sqlite3.Connection:
    conn = sqlite3.connect(db_path)

    conn.execute('''
        CREATE TABLE IF NOT EXISTS gff_data (
            seqID TEXT PRIMARY KEY,
            parent_ID TEXT,
            gene_start INTEGER,
            gene_end INTEGER,
            nuc_length INTEGER,
            aa_length INTEGER,
            strand TEXT,
            COG_ID TEXT,
            cugo_number INTEGER,
            cugo_start TEXT,
            cugo_end TEXT
        )
    ''')

    conn.execute('''
            CREATE TABLE IF NOT EXISTS tmhmm_data (
                seqID TEXT PRIMARY KEY,
                no_tmh INTEGER
            )
        ''')

    conn.execute('CREATE INDEX IF NOT EXISTS idx_gff_seqid ON gff_data(seqID)')
    conn.execute('CREATE INDEX IF NOT EXISTS idx_tmhmm_seqid ON tmhmm_data(seqID)')

    conn.commit()
    return conn


def extract_gene_info(line: list) -> tuple:
    """
    Extracts gene information from a GFF/GTF annotation line.

    Args:
        line (list): parsed annotation line with genomic feature data
        tmhmm_dict (dict): mapping of protein IDs to transmembrane helix counts

    Returns:
        seqID (str): sequence identifier with underscores normalized
        COG_ID (str): COG identifier or 'NA' if not found
        parent (str): chromosome/contig name
        direction (str): strand direction (+/-)
        gene_start (str): start position
        gene_end (str): end position
        nuc_length (int): nucleotide sequence length
        aa_length (int): amino acid sequence length
    """
    # parse annotation field (9th column) for sequence and COG identifiers
    annotation = line[8].split(';')
    seqID = annotation[0].split('=')[1].replace('__', '_')
    COG_ID = annotation[1].split('=')[1] if len(annotation) > 1 else 'NA'

    # extract basic genomic coordinates and features
    parent = line[0]  # chromosome/contig name
    direction = line[6]  # strand direction (+/-)
    gene_start = line[3]  # start position
    gene_end = line[4]  # end position

    # calculate sequence lengths
    nuc_length = abs(int(gene_end) - int(gene_start)) + 1  # nucleotide length
    aa_length = int(nuc_length / 3)  # amino acid length
    return seqID, COG_ID, parent, direction, gene_start, gene_end, nuc_length, aa_length


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
            - COG_ID (str): COG functional category identifier
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
        seqID, COG_ID, parent, direction, gene_start, gene_end, nuc_length, aa_length = extract_gene_info(
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
        reformat_data.append([seqID, parent, gene_start, gene_end, nuc_length, aa_length,
                              direction, COG_ID, cugo, cugo_start, cugo_end])

        # update tracking variables for next iteration
        prev_line = clean_line
        prev_direction = direction
        prev_parent = parent
        prev_cugo = cugo

    # process the last gene in the file (no next gene to compare)
    if len(clean_line) == 9 and clean_line[2] == 'CDS':
        seqID, COG_ID, parent, direction, gene_start, gene_end, nuc_length, aa_length = extract_gene_info(
            clean_line)

        # handle last gene - no next gene, so use none for next_parent
        cugo, cugo_start, cugo_end, cugo_count, cugo_size_count = cugo_boundaries(
            direction, prev_direction, None,  # no next direction for last gene
            parent, prev_parent, None,  # no next parent for last gene
            cugo_count, cugo_size, cugo_size_count, prev_cugo
        )

        # store final gene data (all 11 fields for database)
        reformat_data.append([seqID, parent, gene_start, gene_end, nuc_length, aa_length,
                              direction, COG_ID, cugo, cugo_start, cugo_end])

        # finalize cugo size tracking
        cugo_size[cugo] = cugo_size_count

    # Return complete data for database storage
    return reformat_data


def process_gff_file(filepath: str):
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

        gff_data = parse_single_gff(content)
        return gff_data
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
                    prot_id = parts[0].replace('___', '__')
                    try:
                        no_tmh = int(parts[1])
                    except ValueError:
                        no_tmh = 0
                    tmhmm_data.append((prot_id, no_tmh))
        return tmhmm_data
    except Exception as e:
        logger.warning(f"Failed to process {tmhmm_filepath}: {e}")
        return []


def export_to_tsv(db_path: str, output_tsv: str):
    """Export merged data from SQLite to TSV"""
    conn = sqlite3.connect(db_path)

    # Join GFF and TMHMM data
    query = '''
        SELECT 
            g.seqID,
            g.parent_ID,
            g.aa_length,
            g.strand,
            g.COG_ID,
            g.cugo_number,
            COALESCE(t.no_tmh, 0) as no_TMH
        FROM gff_data g
        LEFT JOIN tmhmm_data t ON g.seqID = t.seqID
        ORDER BY g.seqID
    '''

    cursor = conn.execute(query)

    # Write TSV
    with open(output_tsv, 'w') as f:
        # Header
        f.write('seqID\tparent_ID\taa_length\tstrand\tCOG_ID\tCUGO_number\tno_TMH\n')

        # Data
        for row in cursor:
            f.write('\t'.join(map(str, row)) + '\n')

    conn.close()


def parse(gff_tar_path: str,
          tmhmm_tar_path: str = None,
          output_dir: str = None,
          globdb_version: int = None,
          force: bool = False,
          db_path: str = "temp_genome_data.db",
          cleanup_db: bool = False):
    """
    Main parsing function - single-threaded version for memory efficiency
    """
    output_tsv = ensure_path(output_dir, f'globdb_r{globdb_version}_cugo', force=force)

    logger.info("Processing files single-threaded for memory efficiency")

    if output_dir is None:
        output_dir = '.'

    output_path = Path(output_dir)
    tempdir = output_path / 'tempdir'
    if tempdir.exists():
        shutil.rmtree(tempdir)
    tempdir.mkdir(parents=True, exist_ok=True)

    # Extract tar files
    logger.info("Extracting GFF tar file...")
    subprocess.run(["tar", "-xzf", gff_tar_path, "-C", str(tempdir)], check=True)

    if tmhmm_tar_path:
        logger.info("Extracting TMHMM tar file...")
        subprocess.run(["tar", "-xzf", tmhmm_tar_path, "-C", str(tempdir)], check=True)

    gff_files = list(tempdir.rglob("*_cog.gff.gz"))
    tmhmm_files = list(tempdir.rglob("*_tmhmm_clean"))

    logger.info(f"Found {len(gff_files)} GFF files and {len(tmhmm_files)} TMHMM files")

    # Setup database
    conn = setup_database(db_path)

    # Process GFF files sequentially
    logger.info("Processing GFF files...")
    for gff_file in tqdm(gff_files, desc="GFF files"):
        gff_data = process_gff_file(str(gff_file))
        if gff_data:
            conn.executemany("""
                INSERT OR REPLACE INTO gff_data
                (seqID, parent_ID, gene_start, gene_end, nuc_length,
                 aa_length, strand, COG_ID, cugo_number, cugo_start, cugo_end)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, gff_data)
            conn.commit()

    # Process TMHMM files sequentially
    if tmhmm_files:
        logger.info("Processing TMHMM files...")
        for tmhmm_file in tqdm(tmhmm_files, desc="TMHMM files"):
            tmhmm_data = process_tmhmm_file(str(tmhmm_file))
            if tmhmm_data:
                conn.executemany("""
                    INSERT OR REPLACE INTO tmhmm_data
                    (seqID, no_tmh)
                    VALUES (?, ?)
                """, tmhmm_data)
                conn.commit()

    conn.close()

    # Export to TSV
    logger.info("Exporting to TSV...")
    export_to_tsv(db_path, output_tsv)

    # Cleanup
    if cleanup_db:
        Path(db_path).unlink(missing_ok=True)
    shutil.rmtree(tempdir)

    logger.info(f"Parsing complete. Output saved to: {output_tsv}")

# ======================================
# FUNCTION DEFINITIONS FOR CUGO PLOTTING
# ======================================
def load_cugo_context(path: str):
    return pd.read_csv(path, sep='\t', na_values='', keep_default_na=False)

def extract_flanking_window(df: pd.DataFrame, lower: int, upper: int):
    flank_list = list(range(lower, upper + 1))
    flank_str = [str(i) for i in flank_list]
    return df[flank_str]

def top_context(df: pd.DataFrame, top_n: int):
    """
    Extracts top N most frequent values and their counts from each column.

    Args:
        df (pd.DataFrame): input dataframe to analyze
        top_n (int): number of top values to extract per column

    Returns:
        tuple: (id_df, count_df) where:
            - id_df (pd.DataFrame): top N values for each column
            - count_df (pd.DataFrame): corresponding counts for each value
    """
    top_ids = [[] for _ in range(top_n)]
    top_counts = [[] for _ in range(top_n)]

    for col in df.columns:
        # get value counts for current column
        counts = df[col].value_counts()
        for i in range(top_n):
            if i < len(counts):
                top_ids[i].append(counts.index[i])
                top_counts[i].append(counts.iloc[i])
            else:
                # fill with nan if fewer than top_n values exist
                top_ids[i].append(np.nan)
                top_counts[i].append(np.nan)

    id_df = pd.DataFrame(top_ids, columns=df.columns)
    count_df = pd.DataFrame(top_counts, columns=df.columns)
    return id_df, count_df


def plot_top_cogs_per_position(
        context_path: str,
        flank_lower: int,
        flank_upper: int,
        top_n: int,
        title: str = 'Top COGs per Position',
        save: bool = False,
        ax: plt.Axes = None,
        plot_path: str = None
):
    """
    Creates scatter plot showing top COG categories at each genomic position.

    Args:
        context_path (str): path to context data file
        flank_lower (int): lower bound of flanking region
        flank_upper (int): upper bound of flanking region
        top_n (int): number of top COGs to display per position
        title (str): plot title
        save (bool): whether to save plot to file
        ax (plt.Axes, optional): matplotlib axes object to plot on
        plot_path (str, optional): path to save plot file

    Returns:
        tuple: (positions, boundaries) where:
            - positions (list): list of position values
            - boundaries (list): list of position boundary values
    """
    # load context data and extract cog information
    cont = pd.read_csv(context_path, sep='\t', na_values='', keep_default_na=False)
    extract = cont.loc[cont['COG_ID'] == 'COG_ID']

    # load cog color mapping
    script_dir = Path(__file__).resolve().parent
    color_yaml_path = script_dir / 'cog_colors.yaml'

    with open(color_yaml_path, 'r') as f:
        cog_color_map = yaml.safe_load(f)

    # select columns for specified flanking region
    flank_cols = [str(i) for i in range(flank_lower, flank_upper + 1)]
    cugo_df = extract[flank_cols]
    positions = [int(col) for col in cugo_df.columns]

    # extract top n values and counts for each position
    top_ids = [[] for _ in range(top_n)]
    top_counts = [[] for _ in range(top_n)]

    for col in cugo_df.columns:
        counts = cugo_df[col].value_counts()
        for i in range(top_n):
            if i < len(counts):
                top_ids[i].append(counts.index[i])
                top_counts[i].append(counts.iloc[i])
            else:
                top_ids[i].append(np.nan)
                top_counts[i].append(np.nan)

    id_df = pd.DataFrame(top_ids, columns=cugo_df.columns)
    count_df = pd.DataFrame(top_counts, columns=cugo_df.columns)

    # calculate positioning for multiple ranks at each position
    subtick_width = 0.8 / top_n
    subtick_offset = (top_n - 1) * subtick_width / 2

    x_pos, y_values, cog_labels, point_colors = [], [], [], []
    all_xticks, all_xlabels = [], []

    # create figure if no axes provided
    if ax is None:
        figsize = (max(8, len(positions) * 0.8), 7)
        fig, ax = plt.subplots(figsize=figsize)

    # prepare data points for each position and rank
    for pos_idx, pos in enumerate(positions):
        for rank in range(top_n):
            cog_id = id_df.loc[rank, str(pos)]
            count = count_df.loc[rank, str(pos)]

            # calculate x position with offset for multiple ranks
            x = pos + (rank * subtick_width) - subtick_offset
            x_pos.append(x)
            y_values.append(count)
            all_xticks.append(x)

            if pd.isna(cog_id):
                cog_labels.append('nan')
                point_colors.append('#ffffff')
                all_xlabels.append('nan')
            else:
                cog_labels.append(cog_id)
                point_colors.append(cog_color_map.get(cog_id, '#aaaaaa'))
                all_xlabels.append(cog_id)

    # create scatter plot
    ax.scatter(
        x_pos,
        y_values,
        color=point_colors,
        s=150,
        zorder=3,
        edgecolors='black',
        linewidths=0.7
    )

    # set x-axis labels and ticks
    ax.set_xticks(all_xticks)
    ax.set_xticklabels(all_xlabels, rotation=90, fontsize=14, ha='center')

    # add vertical lines between positions
    for pos in positions[:-1]:
        boundary = pos + 0.5
        ax.axvline(boundary, color='gray', linestyle='--', linewidth=0.7, zorder=1)

    # add position labels at top
    max_y = max([y for y in y_values if not pd.isna(y)], default=0)
    label_offset = max_y * 0.05

    for pos in positions:
        ax.text(pos, max_y + label_offset, str(pos), ha='center', va='bottom',
                fontsize=16, fontweight='bold')

    # set axis limits and labels
    ax.set_ylim(-max_y * 0.05, max_y * 1.15)
    ax.set_xlim(positions[0] - 0.5, positions[-1] + 0.5)

    ax.tick_params(axis='y', labelsize=14)
    ax.set_ylabel('Count', fontsize=16)
    ax.set_title(title, fontsize=18, fontweight='bold')

    # save plot if requested
    if save:
        plt.tight_layout()
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')

    return positions, [pos + 0.5 for pos in positions[:-1]]


def plot_size_per_position(context_path: str,
                           flank_lower: int,
                           flank_upper: int,
                           title: str = 'aa_length Density per Position',
                           save: bool = False,
                           ax: Optional[plt.Axes] = None,
                           pos_boundaries: Optional[list] = None,
                           bin_width: int = 10,
                           y_range: int = None,
                           plot_path: str = None
                           ):
    """
    Creates 1D-density plot showing sequence length distribution across genomic positions.

    Args:
        context_path (str): path to context data file
        flank_lower (int): lower bound of flanking region
        flank_upper (int): upper bound of flanking region
        title (str): plot title
        save (bool): whether to save plot to file
        ax (plt.Axes, optional): matplotlib axes object to plot on
        pos_boundaries (list, optional): position boundaries to draw as vertical lines
        bin_width (int): width of length bins for histogram
        y_range (int, optional): maximum y-axis range to display
        plot_path (str, optional): path to save plot file

    Returns:
        None
    """
    # load context data and extract amino acid length information
    cont = pd.read_csv(context_path, sep='\t', na_values='', keep_default_na=False)
    extract = cont.loc[cont['feat_type'] == 'aa_length']

    # select columns for specified flanking region
    flank_cols = [str(i) for i in range(flank_lower, flank_upper + 1)]
    cugo_df = extract[flank_cols]
    positions = [int(col) for col in cugo_df.columns]

    # determine bin edges for length histograms
    all_lengths = cugo_df.values.flatten()
    all_lengths = all_lengths[~pd.isna(all_lengths)].astype(float)
    max_len = all_lengths.max()
    bin_edges = np.arange(0, max_len + bin_width, bin_width)
    n_bins = len(bin_edges) - 1

    # create histogram data for each position
    heat_data = np.zeros((n_bins, len(flank_cols)))
    position_counts = []

    for col_idx, col in enumerate(flank_cols):
        values = cugo_df[col].dropna().astype(float)
        hist, _ = np.histogram(values, bins=bin_edges)
        heat_data[:, col_idx] = hist
        position_counts.append(len(values))

    # create figure if no axes provided
    if ax is None:
        fig, ax = plt.subplots()

    # set up colormap and normalization
    cmap = get_cmap('Blues')
    norm = Normalize(vmin=0, vmax=heat_data.max())

    # draw heatmap rectangles
    rect_width = 0.8
    for col_idx, pos in enumerate(positions):
        for y in range(heat_data.shape[0]):
            val = heat_data[y, col_idx]
            color = cmap(norm(val))
            ax.add_patch(plt.Rectangle(
                (pos - rect_width / 2, y), rect_width, 1, color=color, linewidth=0
            ))

    # add position boundaries if provided
    if pos_boundaries is not None:
        for boundary in pos_boundaries:
            ax.axvline(boundary, color='gray', linestyle='--', linewidth=0.7, zorder=1)

    # set x-axis labels with position and count information
    ax.set_xlim(positions[0] - 0.5, positions[-1] + 0.5)
    xtick_labels = [f'{pos}\n(n={count})' for pos, count in zip(positions, position_counts)]
    ax.set_xticks(positions)
    ax.set_xticklabels(xtick_labels, fontsize=14)

    # set y-axis range
    if y_range:
        n_bins = int(y_range / bin_width)

    ax.set_ylim(0, n_bins)

    # create y-axis tick labels for length bins
    tick_indices = []
    tick_labels = []

    for i in range(n_bins):
        bin_end = int(bin_edges[i + 1])
        if bin_end % 100 == 0:
            tick_indices.append(i)
            tick_labels.append(f'{bin_end}')

    # fallback tick strategy if no multiples of 100
    if not tick_indices:
        tick_step = max(1, n_bins // 10)
        tick_indices = list(range(0, n_bins, tick_step))
        tick_labels = [f'{int(bin_edges[i + 1])}' for i in tick_indices]

    ax.set_yticks(np.array(tick_indices) + 0.5)
    ax.set_yticklabels(tick_labels, fontsize=16)

    # set axis labels and title
    ax.set_xlabel('Position', fontsize=18)
    ax.set_ylabel(f'length (bin size: {bin_width})', fontsize=18)
    ax.set_title(title, fontsize=20)

    # create colorbar
    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])

    plt.tight_layout()
    if save:
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')

# ======================================
# CUGO COMMAND LINE TOOLS
# ======================================
def context(fasta_path: str,
            cugo_path: str,
            cugo_range: int,
            output_dir: str,
            force: bool = False,
            ):
    """
    Extract genomic context windows around target proteins from CUGO data.

    Args:
        fasta_path: Path to FASTA file
        cugo_path: Path to CUGO file
        cugo_range: Range around target protein for context
        output_dir: Directory for output files
        force: If true, existing files are overwritten

    Returns:
        str: Path to output context file
    """
    # load protein identifiers from fasta or id file
    if fasta_path:
        protein_name = determine_dataset_name(fasta_path, '.', 0)
        output_file = ensure_path(output_dir, f'{protein_name}_context.tsv', force=force)
        sequences = read_fasta_to_dict(fasta_path)
        protein_identifiers = list(sequences.keys())

        # setup special log file
        log_file = ensure_path(output_dir, f'{protein_name}_missing_files.log', force=force)
        logging.basicConfig(
            filename=log_file,
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
    else:
        logger.error('Path to input FASTA file must be provided.')
        raise ValueError('You must provide either a FASTA file or a list of protein IDs.')



    # load cugo data
    with gzip.open(cugo_path, 'rt') as f:
        df = pd.read_csv(f, sep='\t', na_filter=False)

    results = None
    protein_set = set(protein_identifiers)

    # filter for target proteins
    target_df = df[df['seqID'].isin(protein_set)]

    if target_df.empty:
        logging.warning('No target proteins found in the CUGO file.')
        return

    for _, row in target_df.iterrows():
        id = row['seqID']
        target_cugo = int(row['CUGO_number'])
        target_parent = row['parent_ID']
        target_strand = row['strand']

        # extract context window around target protein
        parent_df = df[df['parent_ID'] == target_parent]
        cugo_context = parent_df[
            (parent_df['CUGO_number'] >= (target_cugo - cugo_range)) &
            (parent_df['CUGO_number'] <= (target_cugo + cugo_range))
            ]

        # reverse order for negative strand
        if target_strand == '-':
            cugo_context = cugo_context.iloc[::-1]

        cugo_context = cugo_context.reset_index(drop=True)

        # center index on target protein
        target_index = cugo_context.index[cugo_context.seqID == id].item()
        cugo_context.index = cugo_context.index - target_index
        cugo_context = cugo_context.transpose()

        # combine results
        if results is None:
            results = cugo_context
        else:
            results = pd.concat([results, cugo_context], join='outer', axis=0)

    # save results
    if results is not None:
        results = results.reset_index().rename(columns={'index': 'feat_type'})
        results.to_csv(output_file, sep='\t', index=False)

    return output_file


def cugo_plot(context_path: str,
              flank_lower: int,
              flank_upper: int,
              top_n: int,
              output: str,
              cugo: bool = False,
              size: bool = False,
              all_plots: bool = False,
              bin_width: int = 10,
              y_range: int = None,
              force: bool = False
              ):
    """
    Generate plots for genomic context analysis: COG distribution and protein size.

    Args:
        context_path: Path to context TSV file
        flank_lower: Lower flank boundary for plotting
        flank_upper: Upper flank boundary for plotting
        top_n: Number of top COGs to display
        output: Output directory for plots
        cugo: Whether to generate COG-only plot
        size: Whether to generate size-only plot
        all_plots: Whether to generate combined plot
        bin_width: Bin width for size plots
        y_range: Y-axis range for size plots
        force: Whether to overwrite existing files
    """
    dataset_name = context_path.removesuffix('_context.tsv')

    # generate cog-only plot
    if cugo:
        if not top_n:
            logger.error('top_n is required if cugo is True')
        else:
            cugo_plot_path = ensure_path(output, f'{dataset_name}_cugo_only.png', force=force)
            plot_top_cogs_per_position(context_path=context_path, flank_lower=flank_lower,
                                       flank_upper=flank_upper, top_n=top_n, save=True,
                                       plot_path=cugo_plot_path)

    # generate size-only plot
    if size:
        size_plot_path = ensure_path(output, f'{dataset_name}_size_only.png', force=force)
        plot_size_per_position(context_path=context_path, flank_lower=flank_lower,
                               flank_upper=flank_upper, save=True, bin_width=bin_width,
                               plot_path=size_plot_path)

    # generate combined plot
    if all_plots:
        if not top_n:
            logger.error('top_n is required if cugo is True')
        else:
            all_plot_path = ensure_path(output, f'{dataset_name}_cugo.png')

            # calculate dynamic figure width
            width = max(8, int((flank_upper - flank_lower + 1) * top_n * 0.6))
            figsize = (width, 12)
            fig, (ax1, ax2) = plt.subplots(2, 1,
                                           figsize=figsize,
                                           gridspec_kw={'height_ratios': [3, 2]})

            # plot top cogs per position
            pos_centers, pos_boundaries = plot_top_cogs_per_position(
                context_path=context_path,
                flank_lower=flank_lower,
                flank_upper=flank_upper,
                top_n=top_n,
                title='Top COGs per Position',
                ax=ax1,
            )

            # plot protein size distribution
            plot_size_per_position(
                context_path=context_path,
                flank_lower=flank_lower,
                flank_upper=flank_upper,
                title='aa_length Density per Position',
                ax=ax2,
                pos_boundaries=pos_boundaries,
                bin_width=bin_width,
                y_range=y_range
            )

            ax1.tick_params(labelbottom=True)

            plt.tight_layout()
            plt.subplots_adjust(hspace=0.4)
            plt.savefig(all_plot_path, dpi=300, bbox_inches='tight')

# ======================================
# CUGO COMMAND LINE WORKFLOW
# ======================================

def cugo(cugo_path: str,
         cugo_range: int,
         fasta_path: str,
         output_dir: str,
         flank_lower: int,
         flank_upper: int,
         top_n: int,
         force: bool = False,
         bin_width: int = 10,
         y_range: int = None):
    """
    Complete CUGO workflow: generate context data and create comprehensive plots.

    Args:
        cugo_path: Path to CUGO file
        cugo_range: Range around target protein for context
        output_dir: Directory for output files
        flank_lower: Lower flank boundary for plotting
        flank_upper: Upper flank boundary for plotting
        top_n: Number of top COGs to display
        force: Whether to overwrite existing files
        fasta_path: Optional path to FASTA file
        bin_width: Bin width for size plots
        y_range: Y-axis range for size plots

    Returns:
        tuple: (context_file_path, plot_file_path)
    """
    # generate context data
    context_file = context(
        fasta_path=fasta_path,
        cugo_path=cugo_path,
        cugo_range=cugo_range,
        output_dir=output_dir,
        force=force,
    )

    if context_file is None:
        logger.error("Context generation failed - no plots will be created")
        return None, None

    # create comprehensive plots
    cugo_plot(
        context_path=context_file,
        flank_lower=flank_lower,
        flank_upper=flank_upper,
        top_n=top_n,
        output=output_dir,
        all_plots=True,
        bin_width=bin_width,
        y_range=y_range,
        force=force
    )

    # determine plot file path
    dataset_name = context_file.removesuffix('_context.tsv')
    plot_file = f"{dataset_name}_cugo.png"

    return context_file, plot_file
