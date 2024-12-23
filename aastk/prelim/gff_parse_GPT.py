import pandas as pd
from pathlib import Path


def parse_gff_to_tab(gff_file: str, out_file: str) -> None:
    """
    Parses an Anvi'o-generated GFF file and identifies CUGO (Colocated Unidirectional Gene Organization) information.

    Args:
        gff_file (str): Path to the input GFF file.
        out_file (str): Path to the output tab-delimited file.

    Returns:
        None
    """
    # Validate file paths
    gff_path = Path(gff_file)
    out_path = Path(out_file)

    if not gff_path.is_file():
        raise FileNotFoundError(f"Input GFF file '{gff_file}' not found.")

    # Load and clean GFF content
    column_names = [
        "contig_ID", "source", "feature_type", "start", "stop",
        "score", "strand", "phase", "attributes"
    ]

    # Reading GFF into DataFrame
    gff_data = pd.read_csv(
        gff_path, sep="\t", comment="#", names=column_names, usecols=range(9)
    )

    # Filter CDS entries only
    gff_data = gff_data[gff_data["feature_type"] == "CDS"]

    # Initialize state variables for CUGO detection
    cugo_data = []
    prev_strand = None
    prev_contig = None
    cugo_count = 0
    cugo_size_count = 0

    # Iterate through rows to compute CUGO
    for idx, row in gff_data.iterrows():
        contig = row["contig_ID"]
        strand = row["strand"]

        # Extract attributes
        attributes = dict(attr.split("=", 1) for attr in row["attributes"].split(";") if "=" in attr)
        prot_id = attributes.get("ID", "NA")
        cog_id = attributes.get("COG", "NA")

        if contig != prev_contig or strand != prev_strand:
            # New CUGO detected
            cugo_count += 1
            cugo_size_count = 1
        else:
            # Continuation of the same CUGO
            cugo_size_count += 1

        # Save row data with CUGO info
        cugo_data.append({
            "prot_ID": prot_id,
            "contig_ID": contig,
            "start": row["start"],
            "stop": row["stop"],
            "strand": strand,
            "COG_ID": cog_id,
            "CUGO_number": cugo_count,
            "CUGO_size": cugo_size_count
        })

        # Update previous state
        prev_strand = strand
        prev_contig = contig

    # Create a DataFrame and save output
    cugo_df = pd.DataFrame(cugo_data)
    cugo_df.to_csv(out_path, sep="\t", index=False)

    print(f"CUGO data saved to '{out_file}'.")
