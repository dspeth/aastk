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

### Using conda (recommended)
```bash
git clone https://github.com/dspeth/aastk.git
cd aastk
conda env create -f environment.yml
conda activate aastk
pip install -e .
```

### Using pip
tbd

### From source 
tbd


## Brief descriptions of the tools in AASTK

### Protein alignment score ratio (PASR)
PASR is intended for the creation of comprehensive sequence datasets of homologous proteins, using alignment of all sequences in a query dataset to a (small) seed dataset of sequences of interest using DIAMOND. Thus, the query file is the large dataset containing all proteins to be searched, and the seed is the (preliminary) dataset of homologous proteins to be found in the query data.
The PASR workflow can be run using `aastk pasr`, and the helper tool `pasr_select` can be used to select sequences to include in the final dataset, based on the PASR plot.

More extensive documentation for PASR is available on the [documentation](https://globdb.org/aastk) page and all command line options for PASR are available on the [PASR command line reference page](https://globdb.org/aastk/pasr).     

### Clustering alignment score matrix (CASM)
CASM is designed to investigate the structure of a dataset of homologous proteins, such as a protein superfamily. Many protein (super)families contain members with distinct biochemical or physiological functions. Understanding the structure of a protein (super)family is an essential first step in understanding the functional landscape of a protein (super)family.  CASM clusters sequences by generating an N x n alignment score matrix, by aligning all N sequences in a dataset against a subset of n sequences. T-distributed stochastic neighbourhood embedding (t-SNE) is then used to reduce this matrix to from n to 2 dimensions. Clusters are called using DBSCAN, and the t-SNE results can be annotated with metadata and visualized.
The CASM workflow can be run using `aastk casm`, and the helper tool `pasr_select` can be used to select sequences to include in the final dataset, based on the PASR plot.

More extensive documentation for CASM is available on the [documentation](https://globdb.org/aastk) page and all command line options for CASM are available on the [CASM command line reference page](https://globdb.org/aastk/casm).     

### Colocated unidirectional gene organization (CUGO)
CUGO is intended to retrieve and visualise the consensus genomic context of a dataset of homologous proteins. To do so, it uses information from the AASTK SQL database of the proteins in the GlobDB genomes. Each GlobDB genome is partitioned into CUGO units, that are limited by a strand change of encoded proteins, or when a contig ends. For each sequence in the input file, aastk cugo determines it's CUGO unit, and then extracts all sequences belonging to this CUGO unit, and optionally adjacent CUGOs. These sequences are then used for visualisation.
Consensus genomic context is visualised in three plots combined into one overview figure. The first plot shows the three most common annotations per genomic position, the second the density of amino acid sequence length per position, and the third the density of number of transmembrane helices per position (see figure below). This design ensures that there is no upper limit to the size of the query dataset size from a visualisation perspective, although data retrieval time scales with query size. CUGO uses the AASTK SQL database, which as of release 226 is xxx Gb.

### Metadata retrieval (Meta)
Meta allows for retrieval of sequence metadata from the AASTK SQL database based on a protein fasta file or a list of protein identifiers. There are two types of metadata in the SQL database, those linked to protein sequences directly, and those linked to the genomes encoding the protein sequences. Available metadata categories are annotation (protein linked), taxonomy, culture collection availability, and two levels of environmental metadata (all genome linked). The selected metadata are written to a tsv file.


## Usage examples

### PASR
`aastk pasr -m BLOSUM45 -q query.fasta -s seed.fasta -o output_dir`  
Where:  
`-m` specifies the scoring matrix to be used, can be BLOSUM62 or BLOSUM45  
`-q` specifies the query dataset in fasta format  
`-s` specifies the seed database of homologous proteins, in fasta format  
`-o` specifies the output directory  
All command line options for PASR are available on the [PASR command line reference page](https://globdb.org/aastk/pasr).  

### CASM
`aastk casm --fasta input.faa -o output_directory --subset_size 1000 -n 4 -p 1000`  
Where:  
`--fasta` and `-o` control the input and output   
`--subset_size` determines the number of randomly selected proteins in the subset  
`-n` specifies the number of threads to use  
`-p` is the t-SNE perplexiity  
All command line options for CASM are available on the [CASM command line reference page](https://globdb.org/aastk/casm).  

### CUGO
`aastk cugo -r 0 -l -3 -u 6 -d aastk_sql_database.db -f fasta_file.faa -o cugo_output_dir`  
Where:  
`-r` is the range of CUGO numbers to be considered  
`-l` specifies the number of positions upstream of the gene of interest  
`-u` specifies the number of positions downstream of the gene of interest  
`-d` is the path to the AASTK SQL database  
`-f` and -o control the input and output  
All command line options for CASM are available on the [CUGO command line reference page](https://globdb.org/aastk/cugo).  

### Meta
`aastk meta --all_metadata -d aastk_sql_database.db -f fasta_file.faa -o output_directory`  
Where:  
`--all_metadata` controls which metadata is added to the output tsv file  
`-d` is the path to the AASTK SQL database  
`-f` and `-o` control the input and output  
All command line options for CASM are available on the [Meta command line reference page](https://globdb.org/aastk/meta).  


## Example plots
![test PASR image](https://globdb.org/sites/globdb.org/files/inline-images/mobNAR_bsr.png)
