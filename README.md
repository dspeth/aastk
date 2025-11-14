# Amino acid sequence toolkit

### NOTE: aastk is currently in active development, and not yet ready for broad public use. Please be patient

The amino acid sequence toolkit (aastk) combines a set of tools for analysis of protein amino acid sequences from the the [GlobDB](https://globdb.org/). The GlobDB is the most comprehensive genomic resource of species-representative microbial genomes of Bacteria and Archaea, and contains the largest phylogenetic sequence diversity currently available in a single resource. The aastk is designed to leverage this diversity to create and analyze protein sequence datasets for evolutionary and functional studies.

The aastk is currently comprised of three tools:
- PASR: Database construction using protein alignment score ratio (PASR) analyses.
- CASM: clustering and visuzalization of protein (super)family datasets.
- CUGO: Consensus genomic context visualization



## Installation

### Using conda
```bash
git clone https://github.com/dspeth/aastk.git
cd aastk
conda env create -f environment.yml
conda activate aastk
pip install -e .
```

## Usage examples
PASR

CUGO

CASM

## Usage
```bash
aastk
```
