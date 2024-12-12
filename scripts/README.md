# Genome analysis code

This repository contains some of the genome analysis scripts used for the 1000 Genomes ONT project

## Installation

To create a conda environment with the required tools just use

`git clone https://github.com/1kg-ont-vienna/sv-analysis.git`

`cd sv-analysis/scripts/

`make all`

## Downloading and preparing reference files

The genome files are not included in the repository but can be downloaded using

`cd genome/ && ./prepare_genome.sh`

## Alignment

Alignments for the ONT data to both references are available on the IGSR FTP site.

[https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/)

