# aastk/context.py
import sys
from pathlib import Path
import pandas as pd

def extract_genomic_context(prot_ids, cugo_dir, tmhmm_dir=None, cugo_range=0, out_file=None):
    prot_ids_path = Path(prot_ids)
    cugo_dir_path = Path(cugo_dir)
    tmhmm_dir_path = Path(tmhmm_dir) if tmhmm_dir else None
    out_file_path = Path(out_file)

    # Validate inputs
    if not prot_ids_path.is_file():
        print(f"Error: Protein IDs file '{prot_ids}' does not exist.")
        sys.exit(1)
    if not cugo_dir_path.is_dir():
        print(f"Error: CUGO directory '{cugo_dir}' does not exist.")
        sys.exit(1)
    if tmhmm_dir and not tmhmm_dir_path.is_dir():
        print(f"Error: TMHMM directory '{tmhmm_dir}' does not exist.")
        sys.exit(1)

    cugo_context_all = None
    count = 0

    with open(prot_ids_path, "r") as id_list:
        for line in id_list:
            count += 1
            target_ID = line.strip()
            genome_ID = target_ID.rsplit("_", 1)[0]

            # Load the relevant CUGO file
            cugo_file = cugo_dir_path / f"{genome_ID}_cugo.tab"
            genome_cugo_df = pd.read_csv(cugo_file, sep="\t", na_filter=False)

            if tmhmm_dir:
                # Load the relevant TMHMM file if provided
                tmhmm_file = tmhmm_dir_path / f"{genome_ID}_tmhmm_clean"
                genome_tmhmm_df = pd.read_csv(tmhmm_file, sep="\t", na_filter=False)
                parse_df = pd.merge(genome_cugo_df, genome_tmhmm_df, on="prot_ID", how="left")
            else:
                parse_df = genome_cugo_df

            # Extract genomic context for the target ID
            target_select = parse_df[parse_df["prot_ID"] == target_ID]
            target_cugo = target_select["CUGO_number"].item()
            target_parent = target_select["parent_ID"].item()
            target_strand = target_select["strand"].item()

            parent_df = parse_df[parse_df["parent_ID"] == target_parent]
            cugo_context = parent_df[
                (parent_df["CUGO_number"] >= (target_cugo - cugo_range)) &
                (parent_df["CUGO_number"] <= (target_cugo + cugo_range))
            ]
            if target_strand == "-":
                cugo_context = cugo_context.iloc[::-1]
            cugo_context = cugo_context.reset_index(drop=True)
            target_index = cugo_context.index[cugo_context.prot_ID == target_ID].item()
            cugo_context.index = cugo_context.index - target_index
            cugo_context = cugo_context.transpose()

            # Combine context data across multiple targets
            if cugo_context_all is not None:
                cugo_context_all = pd.concat([cugo_context_all, cugo_context], join="outer", axis=0)
            else:
                cugo_context_all = cugo_context

    cugo_context_all = cugo_context_all.reset_index().rename(columns={"index": "feat_type"})
    cugo_context_all.to_csv(out_file_path, sep="\t", index=False)
    print(f"Successfully extracted genomic context to {out_file}")
