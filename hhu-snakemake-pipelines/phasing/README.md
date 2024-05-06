# Description of the Pipeline

- These are phasing experiments performed using WhatsHap to check the quality of the reads.
- The NYGC study (PMID: 36055201) produced raw genotypes along with the final phased panels. We attempted to phase the raw genotypes using various strategies and consequently compared it against the phased panels to investigate concordance.
- Three phasing strategies were applied:
    - Trio-based phasing: 6 trios are present in the sample set and the genotypes of the trios were used for phasing.
    - Longread-based phasing: For all the samples, phasing was done with only the longread data.
    - Trio Longread phasing: For the 6 trios, the trio data and the longread data was used together to phase.
- The three strategies are then compared to the statistical phasing.

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

The phased VCFs will be created under `results/phased-vcf/<phasing experiment type>/*.vcf`. Along with them, there are files of statistics.

The comparison files are in `results/<trio-comparison | no-trio-comparison>`

Plots can be found under `results/plots/`.

### Reproduction instructions

Under the directory `phasing`, the command `snakemake --use-conda -j <number of cores>` can be used to run the entire pipeline.

Due to the large requirements of certain jobs, a high performance computing architecture will be needed. For running the pipeline in such an architecture, refer to the snakemake documentation (https://snakemake.readthedocs.io/en/v7.18.1/) for further details.