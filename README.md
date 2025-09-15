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
