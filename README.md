# Amino acid sequence toolkit

### NOTE: aastk is currently in active development, and not yet ready for broad public use. Please be patient

The amino acid sequence toolkit (aastk) is a set of tools and approaches to get biological insights from amino acid sequences. Aastk is designed to work with the [GlobDB](https://globdb.org/), a genomic resource comprised of a dereplicated set of species-representative microbial genomes.

Current features of aastk are:
- Protein alignment score ratio (PASR) analyses.
- colocalized unidirection gene organization (CUGO) visualization.

## Installation

### Using conda
```bash
conda env create -f environment.yml
conda activate aastk
pip install -e .
```

## Usage
```bash
aastk
```
