# Amino acid sequence toolkit

## Table of contents
1. [AASTK synopsis](#synopsis)  
2. [Installation](#installation)
3. [Brief descriptions of the tools in AASTK](#tool_description)
4. [Usage examples](#usage_examples)
5. [Example plots](#example_plots)
6. [Citation](#citation)
<br />
<br />

## AASTK synopsis <a name="synopsis"></a>
The amino acid sequence toolkit (AASTK) is a suite of tools that enables construction and analyses of protein sequence datasets from the [GlobDB](https://globdb.org/) genomes. The GlobDB is the most comprehensive genomic resource of species-representative microbial genomes of Bacteria and Archaea, and contains the largest phylogenetic sequence diversity currently available in a single resource. AASTK is designed to leverage the GlobDB to create and analyze protein sequence datasets for evolutionary and functional studies.

AASTK uses a precomputed SQL database that includes the protein sequences, their genomic context, functional annotations from KEGG, COG, and PFAM, and metadata on the genomes the protein sequences. This SQL database is available for download on the [GlobDB website](https://globdb.org/downloads).

AASTK currently consists of four tools:  
`pasr` - Creating and updating comprehensive datasets of proteins of interest.  
`casm` - Clustering datasets of complex functionally diverse protein superfamilies.  
`cugo` - Assessing consensus genomic context of a protein dataset.  
`meta` - Retrieving metadata such as taxonomy, environment, and annotations from the AASTK SQL database

Installation instructions, usage examples, and short descriptions of the tools in AASTK can be found below. A more extended description can be found in the [documentation](https://globdb.org/aastk). A full description of all command line tools and arguments is available on command line reference pages for [PASR](https://globdb.org/aastk/pasr), [CASM](https://globdb.org/aastk/casm), [CUGO](https://globdb.org/aastk/cugo), and [Meta](https://globdb.org/aastk/meta), also accessible via the dropdown menu on the AASTK [documentation](https://globdb.org/aastk) page.
<br />
<br />
<br />
  
## Installation <a name="installation"></a>

### Using conda (recommended)
We recommend installing AASTK using conda, in a dedicated environment created for the software.
```bash
conda create aastk
conda activate aastk
conda install -c bioconda aastk
```

In addition to the software and dependencies, AASTK requires a SQL database for full functionality, and is best used with the GlobDB protein complement. A fastA file of the GlobDB proteins can be exported from the AASTK SQL database using `aastk export_fasta` so only a single download is required. 
To set up the AASTK SQL database and GlobDB protein dataset, run the following commands:
```bash
wget https://fileshare.lisc.univie.ac.at/globdb/globdb_r226/globdb_r226_aastk.db.gz
gunzip globdb_r226_aastk.db.gz
aastk export_fasta -d globdb_r226_aastk.db -n 4 -o protein_fasta_dir
```
<br />

### Using pip
```bash
pip install aastk
```
Dependencies:
 - [Python](https://www.python.org/) >= 3.13
 - [DIAMOND](https://github.com/bbuchfink/diamond)
 - [SeqKit](https://github.com/shenwei356/seqkit)

DIAMOND and SeqKit need to be available in your PATH.  
After installation of the AASTK software and dependencies, the AASTK SQL database and GlobDB protein fastA file can be set up as described under 'Using conda' above.
<br />

### From source
```bash
wget https://github.com/dspeth/aastk/archive/refs/tags/v0.1.0.tar.gz
tar -xzf v0.1.0.tar.gz
cd aastk-0.1.0
pip install .
```

Dependencies:
 - [Python](https://www.python.org/) >= 3.13
 - [DIAMOND](https://github.com/bbuchfink/diamond)
 - [SeqKit](https://github.com/shenwei356/seqkit)

DIAMOND and SeqKit need to be available in your PATH.  
After installation of the AASTK software and dependencies, the AASTK SQL database and GlobDB protein fastA file can be set up as described under 'Using conda' above.
<br />
<br />
<br />

## Brief descriptions of the tools in AASTK <a name="tool_description"></a>

### Protein alignment score ratio (PASR)
PASR is intended for the creation of comprehensive sequence datasets of homologous proteins, using alignment of all sequences in a query dataset to a (small) seed dataset of sequences of interest using DIAMOND. Thus, the query file is the large dataset containing all proteins to be searched, and the seed is the (preliminary) dataset of homologous proteins to be found in the query data.
The PASR workflow can be run using `aastk pasr`, and the helper tool `pasr_select` can be used to select sequences to include in the final dataset, based on the PASR plot.

More extensive documentation for PASR is available on the [documentation](https://globdb.org/aastk) page and all command line options for PASR are available on the [PASR command line reference page](https://globdb.org/aastk/pasr).     
<br />

### Clustering alignment score matrix (CASM)
CASM is designed to investigate the structure of a dataset of homologous proteins, such as a protein superfamily. Many protein (super)families contain members with distinct biochemical or physiological functions. Understanding the structure of a protein (super)family is an essential first step in understanding the functional landscape of a protein (super)family.  CASM clusters sequences by generating an N x n alignment score matrix, by aligning all N sequences in a dataset against a subset of n sequences. T-distributed stochastic neighbourhood embedding (t-SNE) is then used to reduce this matrix to from n to 2 dimensions. Clusters are called using DBSCAN, and the t-SNE results can be annotated with metadata and visualized.
The CASM workflow can be run using `aastk casm`, and the helper tool `pasr_select` can be used to select sequences to include in the final dataset, based on the PASR plot.

More extensive documentation for CASM is available on the [documentation](https://globdb.org/aastk) page and all command line options for CASM are available on the [CASM command line reference page](https://globdb.org/aastk/casm).     
<br />

### Colocated unidirectional gene organization (CUGO)
CUGO is intended to retrieve and visualise the consensus genomic context of a dataset of homologous proteins. To do so, it uses information from the AASTK SQL database of the proteins in the GlobDB genomes. Each GlobDB genome is partitioned into CUGO units, that are limited by a strand change of encoded proteins, or when a contig ends. For each sequence in the input file, aastk cugo determines it's CUGO unit, and then extracts all sequences belonging to this CUGO unit, and optionally adjacent CUGOs. These sequences are then used for visualisation.
Consensus genomic context is visualised in three plots combined into one overview figure. The first plot shows the three most common annotations per genomic position, the second the density of amino acid sequence length per position, and the third the density of number of transmembrane helices per position (see figure below). This design ensures that there is no upper limit to the size of the query dataset size from a visualisation perspective, although data retrieval time scales with query size. CUGO uses the AASTK SQL database, which as of release 226 is 310 Gb.
<br />

### Metadata retrieval (Meta)
Meta allows for retrieval of sequence metadata from the AASTK SQL database based on a protein fasta file or a list of protein identifiers. There are two types of metadata in the SQL database, those linked to protein sequences directly, and those linked to the genomes encoding the protein sequences. Available metadata categories are annotation (protein linked), taxonomy, culture collection availability, and two levels of environmental metadata (all genome linked). The selected metadata are written to a tsv file.
<br />
<br />
<br />

## Usage examples <a name="usage_examples"></a>

### PASR
`aastk pasr -m BLOSUM45 -q query.fasta -s seed.fasta -o output_dir`  
Where:  
`-m` specifies the scoring matrix to be used, can be BLOSUM62 or BLOSUM45  
`-q` specifies the query dataset in fasta format  
`-s` specifies the seed database of homologous proteins, in fasta format  
`-o` specifies the output directory  
All command line options for PASR are available on the [PASR command line reference page](https://globdb.org/aastk/pasr).  
<br />

### CASM
`aastk casm --fasta input.faa -o output_directory --subset_size 1000 -n 4 -p 1000`  
Where:  
`--fasta` and `-o` control the input and output   
`--subset_size` determines the number of randomly selected proteins in the subset  
`-n` specifies the number of threads to use  
`-p` is the t-SNE perplexiity  
All command line options for CASM are available on the [CASM command line reference page](https://globdb.org/aastk/casm).  
<br />

### CUGO
`aastk cugo -r 0 -l -3 -u 6 -d aastk_sql_database.db -f fasta_file.faa -o cugo_output_dir`  
Where:  
`-r` is the range of CUGO numbers to be considered  
`-l` specifies the number of positions upstream of the gene of interest  
`-u` specifies the number of positions downstream of the gene of interest  
`-d` is the path to the AASTK SQL database  
`-f` and -o control the input and output  
All command line options for CASM are available on the [CUGO command line reference page](https://globdb.org/aastk/cugo).  
<br />

### Meta
`aastk meta --all_metadata -d aastk_sql_database.db -f fasta_file.faa -o output_directory`  
Where:  
`--all_metadata` controls which metadata is added to the output tsv file  
`-d` is the path to the AASTK SQL database  
`-f` and `-o` control the input and output  
All command line options for CASM are available on the [Meta command line reference page](https://globdb.org/aastk/meta).  
<br />
<br />
<br />

## Example plots <a name="example_plots"></a>
Examples of the graphical output of PASR, CUGO, and CASM. For a full overview of the output files generated by these tools, consult the [AASTK documentation page](https://globdb.org/aastk).
<br />
  
#### Example PASR plot
![Example PASR plot, described in the legend below](https://globdb.org/sites/globdb.org/files/inline-images/mobNAR_bsr.png)
PASR plot of the catalytic, molybdenum containing subunit of an enzyme in the MopB superfamily. Each dot indicates a protein sequence, with the x-axis representing the calculated maximum alignment score, and the y-axis the alignment score against the seed dataset. Color of the dots represents sequence identity of the best hit in the seed dataset. Dots on the 1:1 line represent sequences already present in the seed dataset. Full length sequences (for this dataset) have a calculated maximum score (x-axis) of approx. 4500-5000. Dots with lower values on the x-axis represent partial sequences, either pseudogenes or (more commonly) sequences encoded at the edge of contigs in fragmented genomes. The x-axis is cut off at 150% of the maximum value of the y-axis
<br />
  
#### Example CUGO plot
![Example CUGO plot, described in the legend below](https://globdb.org/sites/globdb.org/files/inline-images/globdb_r226_narG_COG_ID_cugo.png)
The CUGO visualisation consists of three plots. The first plot shows the three most prevalent annotations per position, using a histogram like graph. Colors for each COG are consistent within and across plots. The second plot shows the number of proteins in length bins (default 50 amino acids) position, showing whether proteins at a position are conserved in an annotation independent way. The third plot shows the predicted transmembrane helices at each position. 
<br />
  
#### Example CASM plot after early exaggeration phase.
![Example CASM plot after early exaggeration phase, described in the legend below](https://globdb.org/sites/globdb.org/files/inline-images/GTDB_r202_plus_GEMOTU_molyb_align_tsne_early_clusters.png)
tSNE plot showing 166,445 sequences of the mopB superfamily after early exaggeration, when the clustering with DBSCAN is done. Points are colored by cluster affiliation. 
<br />
  
#### Example of final CASM plot
![Example of the final CASM plot, described in the legend below](https://globdb.org/sites/globdb.org/files/inline-images/GTDB_r202_plus_GEMOTU_molyb_align_tsne_final_clusters.png)
Final tSNE plot showing 166,445 sequences of the mopB superfamily. Points are colored by cluster affiliation, but can also be colored by information from the AASTK SQL database, including taxonomy, environment or culture availability.
<br />
<br />
<br />

## Citation <a name="citation"></a>
There is no publication describing AASTK yet, so please cite this repository when you use AASTK.  
  
In addition, several parts of the software were developed independently and should be credited.

- If you use AASTK with the GlobDB protein dataset, please cite:  
Speth et al. (2025) __GlobDB: a comprehensive species-dereplicated microbial genome resource__  
https://doi.org/10.1093/bioadv/vbaf280

- If you use `aastk pasr`, please cite:  
Speth and Orphan (2018) __Metabolic marker gene mining provides insight in global mcrA diversity and, coupled with targeted genome reconstruction, sheds further light on metabolic potential of the Methanomassiliicoccales__  
https://doi.org/10.7717/peerj.5614

- The environmental data from `aastk meta` is derived from the MetaCoOc software. A manuscript is in preparation, but in the meantime please cite:   
https://github.com/bcoltman/metacooc
