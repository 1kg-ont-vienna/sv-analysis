# Description of the Pipeline

- The long read DNA sequences of the sample cohort can be tagged as "H1" or "H2" using `whatshap`.
- The haplotags were determined using the NYGC phased panels (PMID: 36055201) where 967 samples intersect between their study and our study.
- Due to WhatsHap's limitation at the pseudo-autosomal regions of the X and Y chromosomes, certain post processing steps were added to accurately tag the reads.

## System Requirements

- python=3.7.12
- pandas=1.3.5
- snakemake=7.18.1
- conda=22.9.0

## Installation Guide

### Instructions

The pipeline is run with snakemake and the dependencies of each rule is resolved by the specified conda environment. The pipeline does not require any installation besides installing the overall dependencies specified under *System Requirements*

### Typical install time

Very fast since it only requires installing the dependencies stated above.

## Instructions for use

### Input

Running the pipeline on our data will require downloading the entire dataset. Refer to the *config.yaml* file where it states the required data for the pipeline and where to download them. The required directories and files need to be specified in the variables and then the pipeline can be run.

### Output

The pipeline outputs a set of haplotags as a TSV file under `results/GRCh38/<sample>/<sample>.tsv`.

### Reproduction instructions

Under the directory `haplotagging`, the command `snakemake --use-conda -j <number of cores>` can be used to run the entire pipeline.

Due to the large requirements of certain jobs, a high performance computing architecture will be needed. For running the pipeline in such an architecture, refer to the snakemake documentation (https://snakemake.readthedocs.io/en/v7.18.1/) for further details.