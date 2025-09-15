# Overview

This script automates sequence-homology-based gene prediction via miniprot. See Li et al. 2023 for more details (https://doi.org/10.1093/bioinformatics/btad014).

The pipeline runs through the following main steps:
1. Convert .fasta genome assemblies into miniprot index format.
2. Run miniprot to obtain gene predictions in .gff format
3. Extract nucleotide sequences of predicted genes using AGAT, a .gff file processor
4. Sort (or "combine") individual coding sequences into fasta files by gene (suitable for downstream programs such as PAML).

An additional command separate from the main pipeline filters out poor-quality predictions from the combined fasta files. 

# Requirements

The repository can be cloned with the following command:
```git clone https://github.com/tonyxieuoft/miniprot_pipeline.git```

The python package `Bio` must be installed (details can be found here: https://biopython.org/wiki/Download). 

Other installation requirements include:
- miniprot 
- AGAT, a .gff file processor

We recommend downloading the above as packages in a conda environment (see https://anaconda.org/bioconda/miniprot, https://anaconda.org/bioconda/agat for more details). 

# Usage

To run the main pipeline, `cd` into the `miniprot_pipeline` folder and enter the following command:
```
python3 run_miniprot_pipeline -f /path/to/fasta/genome/assemblies -r /path/to/reference/sequence/directory -d /output/path [OPTIONAL] -i /path/to/miniprot/genome/indexes
```

Each input directory has specific formatting requirements. 

The fasta genome assembly directory specified by the `-f` flag requires the following:
-At least one child directory.
-Child directories contain .fasta genome assemblies named by species
-Each child directory is named after the taxon of interest encompassing the species in the directory.

The reference sequence directory specified by the `-r` flag requires the following: 
-At least one child directory.
-Child directories contain reference sequence fasta files (nucleotide format) named by gene or variant.
-Each child directory is named after a taxon of interest in the fasta genome assembly directory. The reference sequences will be used query the genome assemblies in the directory with the matching name. 

As indicated, the `-i` is optional and requires the prior creation of miniprot indexes. It is only meant to speed up the pipeline in subsequent instances after already running it once.

Sample input directories are provided for reference. 
