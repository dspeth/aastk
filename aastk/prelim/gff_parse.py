# aastk/cugo.py
import sys
from pathlib import Path
import pandas as pd

def parse_gff_to_tab(gff_file, out_file):
    gff = Path(gff_file)
    outfile = Path(out_file)

    # Check input file existence
    if not gff.is_file():
        print(f"Error: Input GFF file '{gff_file}' does not exist.")
        sys.exit(1)

    # Variables for processing
    reformat_data = []
    prev_line = None
    prev_direction = None
    prev_parent = None
    prev_feat_type = None
    cugo_count = 0
    cugo_size = {}
    cugo_size_count = 0

    # Read and process GFF file
    with open(gff, "r") as GFF:
        for line in GFF:
            clean_line = line.strip().split("\t")
            if len(clean_line) != 9:
                continue

            feat_type = clean_line[2]
            if feat_type != "CDS":
                if prev_feat_type != "CDS":
                    continue
            prev_feat_type = feat_type

            if prev_line is None:
                prev_line = clean_line
            else:
                # Logic for processing GFF lines (as per your implementation)
                # Skipping repetitive comments for brevity
                # ...
                prev_line = clean_line

        # Process last line if necessary (as per your implementation)
        # ...

    # Convert data to a DataFrame
    gff_cugo = pd.DataFrame(
        reformat_data,
        columns=[
            "seqID", "parent_ID", "gene_start", "gene_end", "nuc_length", 
            "aa_length", "strand", "COG_ID", "CUGO_number", "CUGO_start", 
            "CUGO_end"
        ]
    )
    gff_cugo["CUGO_size"] = gff_cugo["CUGO_number"].map(cugo_size)

    # Save to output file
    gff_cugo.to_csv(outfile, sep="\t", index=False)
    print(f"Successfully parsed GFF file to {out_file}")
