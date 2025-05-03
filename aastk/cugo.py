#!/usr/bin/env python3

from .util import ensure_path

import pandas as pd

def extract_gene_info(line):
    """
    Extracts gene info from a GFF line
    Args:
        line (list): A tab-separated line from GFF file split into fields

    Returns:
        tuple: (seqID, COG_ID, parent, direction, gene_start, gene_end, nuc_length, aa_length)
        seqID (str): Sequence identifier, cleaned by replacing ___ with _
        COG_ID (str): COG identifier or "NA" if not present
        parent (str): Parent sequence/contig identifier
        direction (str): Strand direction ("+" or "-")
        gene_start (str): Start position of the gene
        gene_end (str): End position of the gene
        nuc_length (int): Nucleotide length of the gene
        aa_length (int): Amino acid length (nucleotide length / 3)
    """
    annotation = line[8].split(';')
    seqID = annotation[0].split('=')[1].replace('__', '_')
    COG_ID = annotation[1].split('=')[1] if len(annotation) > 1 else 'NA'
    parent = line[0]
    direction = line[6]
    gene_start = line[3]
    gene_end = line[4]
    nuc_length = abs(int(gene_end) - int(gene_start)) + 1
    aa_length = int(nuc_length / 3)
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
    cugo_start, cugo_end = "NA", "NA"

    # middle of CUGO
    if (direction == prev_direction == next_direction) and (parent == prev_parent == next_parent):
        cugo = prev_cugo
        cugo_size_count += 1

    # start of contig/scaffold
    elif parent != prev_parent:
        cugo = cugo_count
        cugo_size_count = 1

        if direction == "+":
            cugo_start = "sequence_edge"
            if parent != next_parent:
                cugo_end = "sequence_edge"
                cugo_size[cugo] = cugo_size_count
                cugo_size_count = 0
                cugo_count += 1
            elif direction != next_direction:
                cugo_end = "strand_change"
                cugo_size[cugo] = cugo_size_count
                cugo_size_count = 0
                cugo_count += 1
        else:  # direction == "-"
            cugo_end = "sequence_edge"
            if parent != next_parent:
                cugo_start = "sequence_edge"
                cugo_size[cugo] = cugo_size_count
                cugo_size_count = 0
                cugo_count += 1
            elif direction != next_direction:
                cugo_start = "strand_change"
                cugo_size[cugo] = cugo_size_count
                cugo_size_count = 0
                cugo_count += 1

    # end of contig/scaffold
    elif parent != next_parent:
        if direction == prev_direction:
            cugo = prev_cugo
            cugo_size_count += 1
            if direction == "+":
                cugo_end = "sequence_edge"
            else:
                cugo_start = "sequence_edge"
            cugo_size[cugo] = cugo_size_count
            cugo_size_count = 0
            cugo_count += 1
        else:
            cugo = cugo_count
            cugo_size_count = 1
            if direction == "+":
                cugo_start = "strand_change"
                cugo_end = "sequence_edge"
            else:
                cugo_start = "sequence_edge"
                cugo_end = "strand_change"
            cugo_size[cugo] = cugo_size_count
            cugo_size_count = 0
            cugo_count += 1

    # Strand change
    elif direction != prev_direction or direction != next_direction:
        if direction != prev_direction and direction == next_direction:
            cugo = cugo_count
            cugo_size_count = 1
            if direction == "+":
                cugo_start = "strand_change"
            else:
                cugo_end = "strand_change"
        elif direction == prev_direction and direction != next_direction:
            cugo = prev_cugo
            cugo_size_count += 1
            if direction == "+":
                cugo_end = "strand_change"
            else:
                cugo_start = "strand_change"
            cugo_size[cugo] = cugo_size_count
            cugo_size_count = 0
            cugo_count += 1
        else:  # direction != prev_direction and direction != next_direction
            cugo = cugo_count
            cugo_size_count = 1
            cugo_start = cugo_end = "strand_change"
            cugo_size[cugo] = cugo_size_count
            cugo_size_count = 0
            cugo_count += 1

    return cugo, cugo_start, cugo_end, cugo_count, cugo_size_count

def process_last_gene(line,
                      prev_direction: str,
                      prev_parent: str,
                      cugo_count: int,
                      cugo_size_count: int,
                      prev_cugo: int,
                      cugo_size: dict):
    seqID, COG_ID, parent, direction, gene_start, gene_end, nuc_length, aa_length = extract_gene_info(line)
    cugo_start, cugo_end = "NA", "NA"

    if (direction == prev_direction) and (parent == prev_parent):
        cugo = prev_cugo
        cugo_size_count += 1
        cugo_size[cugo] = cugo_size_count
        if direction == "+":
            cugo_end = "sequence_edge"
        else:
            cugo_start = "sequence_edge"
    elif parent != prev_parent:
        cugo = cugo_count
        cugo_size_count = 1
        cugo_size[cugo] = cugo_size_count
        cugo_start = cugo_end = "sequence_edge"
    elif direction != prev_direction:
        cugo = cugo_count
        cugo_size_count = 1
        cugo_size[cugo] = cugo_size_count
        if direction == "+":
            cugo_start, cugo_end = "strand_change", "sequence_edge"
        else:
            cugo_start, cugo_end = "sequence_edge", "strand_change"

    return [seqID, parent, gene_start, gene_end, nuc_length, aa_length,
            direction, COG_ID, cugo, cugo_start, cugo_end], cugo_size

def parse(gff_file_path: str,
          output_dir: str):
    reformat_data = []
    cugo_size = {}
    cugo_count = prev_direction = prev_parent = prev_feat_type = prev_cugo = 0
    cugo_size_count = 0
    prev_line = None
    line_count = 0
    cds_count = 0

    # name eventual output file
    identifier = gff_file_path.split('.')[0] + '_cugo'
    output_path = ensure_path(output_dir, identifier)

    try:
        with open(gff_file_path, "r") as GFF:
            for line_count, line in enumerate(GFF, 1):
                clean_line = line.strip().split("\t")

                # skip invalid lines and non-CDS features
                if len(clean_line) != 9:
                    continue

                feat_type = clean_line[2]
                if feat_type == "CDS":
                    cds_count += 1

                if feat_type != "CDS" and prev_feat_type != "CDS":
                    continue
                prev_feat_type = feat_type

                # initialize with first valid line
                if prev_line is None:
                    prev_line = clean_line
                    continue

                # extract information from current line
                seqID, COG_ID, parent, direction, gene_start, gene_end, nuc_length, aa_length = extract_gene_info(
                    prev_line)

                # get context
                next_parent = clean_line[0]
                next_direction = clean_line[6]

                # Determine CUGO membership
                cugo, cugo_start, cugo_end, cugo_count, cugo_size_count = cugo_boundaries(
                    direction, prev_direction, next_direction,
                    parent, prev_parent, next_parent,
                    cugo_count, cugo_size, cugo_size_count, prev_cugo
                )

                # Add to results
                reformat_data.append([seqID, parent, gene_start, gene_end, nuc_length, aa_length,
                                      direction, COG_ID, cugo, cugo_start, cugo_end])

                # Update for next iteration
                prev_line = clean_line
                prev_direction = direction
                prev_parent = parent
                prev_cugo = cugo

            # Process last line if it's valid
            if len(clean_line) == 9 and clean_line[2] == "CDS":
                last_line, cugo_size = process_last_gene(
                    clean_line, prev_direction, prev_parent,
                    cugo_count, cugo_size_count, prev_cugo, cugo_size
                )
                reformat_data.append(last_line)

        # Create dataframe
        gff_cugo = pd.DataFrame(
            reformat_data,
            columns=["seqID", "parent_ID", "gene_start", "gene_end", "nuc_length",
                     "aa_length", "strand", "COG_ID", "CUGO_number", "CUGO_start", "CUGO_end"]
        )
        gff_cugo["CUGO_size"] = gff_cugo["CUGO_number"].map(cugo_size)

        gff_cugo.to_csv(output_path, sep="\t", index=False)

        return gff_cugo

    except Exception as e:
        raise Exception(f"Error processing GFF file: {str(e)}")






