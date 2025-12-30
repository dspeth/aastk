# Amino acid sequence toolkit

### NOTE: aastk is currently in active development, and not yet ready for broad public use. Please be patient
## AASTK synopsis
The amino acid sequence toolkit (AASTK) is a suite of tools that enables construction and analyses of protein sequence datasets from the [GlobDB](https://globdb.org/) genomes. The GlobDB is the most comprehensive genomic resource of species-representative microbial genomes of Bacteria and Archaea, and contains the largest phylogenetic sequence diversity currently available in a single resource. AASTK is designed to leverage the GlobDB to create and analyze protein sequence datasets for evolutionary and functional studies.

AASTK uses a precomputed SQL database that includes the protein sequences, their genomic context, functional annotations from KEGG, COG, and PFAM, and metadata on the genomes the protein sequences. This SQL database is available for download on the [GlobDB website](https://globdb.org/downloads).

AASTK currently consists of four tools:  
`pasr` - Creating and updating comprehensive datasets of proteins of interest.  
`casm` - Clustering datasets of complex functionally diverse protein superfamilies.  
`cugo` - Assessing consensus genomic context of a protein dataset.  
`meta` - Retrieving metadata such as taxonomy, environment, and annotations from the AASTK SQL database

Installation instructions, usage examples, and short descriptions of the tools in AASTK can be found below. A more extended description can be found in the [documentation](https://globdb.org/aastk). A full description of all command line tools and arguments is available on command line reference pages for [PASR](https://globdb.org/aastk/pasr), [CASM](https://globdb.org/aastk/casm), [CUGO](https://globdb.org/aastk/cugo), and [Meta](https://globdb.org/aastk/meta), also accessible via the dropdown menu on the AASTK [documentation](https://globdb.org/aastk) page.


## Installation

#### Using conda (recommended)
```bash
git clone https://github.com/dspeth/aastk.git
cd aastk
conda env create -f environment.yml
conda activate aastk
pip install -e .
```

#### Using pip
tbd

#### From source 
tbd


## Usage examples

#### Protein alignment score ratio (PASR)
PASR is intended for the creation of comprehensive sequence datasets of homologous proteins, using alignment of all sequences in a query dataset to a (small) seed dataset of sequences of interest using DIAMOND. The PASR workflow can be run using `aastk pasr`, as shown in the usage examples below:  
  
All command line options for PASR are available on the [PASR command line reference page](https://globdb.org/aastk/pasr).     
<br />
#### Clustering alignment score matrix (CASM)
CASM is designed to investigate the structure of a dataset of homologous proteins, such as a protein superfamily. Many protein (super)families contain members with distinct biochemical or physiological functions. Understanding the structure of a protein (super)family is an essential first step in understanding the functional landscape of a protein (super)family. Existing approaches such as phylogenetic analyses or sequence similarity networks don't scale well to the size of modern datasets.  
CASM clusters sequences by generating an N x n alignment score matrix, by aligning all N sequences in a dataset against a subset of n sequences. T-distributed stochastic neighbourhood embedding (t-SNE) is then used to reduce this matrix to from n to 2 dimensions. Clusters are called using DBSCAN, and the t-SNE results can be annotated with metadata and visualized.
The CASM workflow can be run using `aastk casm`, as shown in the usage examples below:  


All command line options for CASM are available on the [CASM command line reference page](https://globdb.org/aastk/casm).     
<br />
####Colocated unidirectional gene organization (CUGO)
