from .util import *

import sys
from pathlib import Path
import random
import numpy as np
import subprocess
import pandas as pd

def fasta_subsample(fasta: str, output_path: str, subset_size: int):
    sequences = read_fasta_to_dict(fasta)
    subset_keys = random.sample(sequences.keys(), subset_size)
    subset_dict = {k: sequences[k] for k in subset_keys}
    with open(output_path, 'w') as f:
        for header, seq in sequences.items():
            f.write(f">{header}\n")
            f.write(f'{seq}\n')



def asm_clust(fasta: str, output: str, subset_size: int):
    fasta_subsample(fasta, output, subset_size)
