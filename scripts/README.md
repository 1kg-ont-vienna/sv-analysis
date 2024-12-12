# Genome analysis code

This repository contains some of the genome analysis scripts used for the 1000 Genomes ONT project

## Installation

To create a mamba environment with the required tools just use

`git clone https://github.com/1kg-ont-vienna/sv-analysis.git`

`cd sv-analysis/scripts/`

`make all`

## Downloading and preparing reference files

The genome files are not included in the repository but can be downloaded using

`cd genome/ && ./prepare_genome.sh`

## Alignment

Alignments for the ONT data to both references are available on the IGSR FTP site.

[https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/)

The scripts assume these alignments are available in the subfolders `./hg38/` and `./t2t/`.

## Data release

You will also need a copy of the latest data release for the 1000 Genomes ONT project to run some of the scripts.

[https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/release/](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/release/)

## Alignment statistics

To calculate the alignment error rate, genome coverage and other statistics.

`cd qc/ && ./qc.sh`

## Lifting of SVs to GRCh38

This workflow uses the phased SVs to implant SVs into CHM13 haplotypes, which are then aligned to GRCh38 and called via [svim-asm](https://github.com/eldariont/svim-asm).

