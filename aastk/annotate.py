"""
annotate.py – Step 1: Gene calling with Pyrodigal
Converts a genome FASTA into predicted amino acid sequences (.faa)
"""

import pyrodigal
from Bio import SeqIO

# --- Input and output (set manually) ---
genome_fasta = "test_data/GCF_000005845.2_ASM584v2_genomic.fna"
output_fasta = "predicted_proteins.faa"

# --- Read genome sequence ---
records = list(SeqIO.parse(genome_fasta, "fasta"))

# --- Initialize gene finder ---
orf_finder = pyrodigal.GeneFinder(meta=True)



# --- predict genes and translation---
with open(output_fasta, "w") as out_faa:
    for record in records:
        results = orf_finder.find_genes(str(record.seq))
        for i , gene in enumerate (results) :
            header = f">{record.id}_gene{i+1} {gene.begin}:{gene.end} ({'+' if gene.strand == 1 else '-'})"
            aa_seq = gene.translate()
            out_faa.write(f"{header} \n {aa_seq}\n")

print (f"Predicted proteins written to {output_fasta}")