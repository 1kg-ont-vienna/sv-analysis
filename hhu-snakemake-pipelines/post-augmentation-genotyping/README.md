# Description of the Pipeline

- Uses the pseudohaplotypes from the SAGA pipeline, HPRC assemblies, and the rGFA file (which has the pseudohaplotypes included) to create genotype calls from Giggles.
- Creates the panel VCF by aligning assemblies (and pseudohaplotypes) back to the rGFA and finding the alleles from each bubble.
- Aligns the long read DNA sequences to the rGFA using minigraph and then sorts it with gaftools.
- Genotypes the samples from the callset using the panel VCF as the set of genotypable variant positions using Giggles.
- Post-processes some basic graphs and statistics from the genotypes (most of the post-genotyping analysis is shown in the folder `hhu-analysis-scripts`)

## System Requirements

- python=3.7.12
- pandas=1.3.5
- snakemake=7.18.1
- conda=22.9.0
- bgzip=1.9
- tabix=1.9
- minigraph=0.20-r559
- gaftools (commit ID: [903921ccc7ef4ffb6dd150974c080e75b5d72aac](https://github.com/marschall-lab/gaftools/tree/903921ccc7ef4ffb6dd150974c080e75b5d72aac))
- giggles=1.0

## Installation Guide

### Instructions

The pipeline is run with snakemake and the dependencies of each rule is resolved by the specified conda environment.

The following packages need installation:

- minigraph
    - link to github: https://github.com/lh3/minigraph
    - the path to the executable needs to be provided

- gaftools
    - link to github: https://github.com/marschall-lab/gaftools
    - for this pipeline, gaftools needs to be installed in a conda environment called `gaftools-env`. Instructions on conda installation is provided in the github page.

- giggles
    - link to github: https://github.com/samarendra-pani/giggles
    - for this pipeline, giggles needs to be installed in a conda environment called `giggles-env`. Instructions on conda installation is provided in the github wiki.

## Instructions for use

### Input

Running the pipeline on our data will require downloading the entire dataset. Refer to the *config.yaml* file where it states the required data for the pipeline and where to download them. The required directories and files need to be specified in the variables and then the pipeline can be run.

### Output

The main outputs of this pipeline are:

- sorted GAF reads at `results/<callset>/ont-alignments/<sample>.sorted.gaf.gz`
- the VCF panel used as input for giggles at `results/<callset>/panel/giggles-ready_multiallelic.vcf.gz` and accompanying VCFs.
- giggles genotypes in multiallelic (`results/<callset>/genotypes/multisample-multiallelic.vcf.gz`) and biallelic form (`results/<callset>/genotypes/multisample-biallelic.vcf.gz`)

### Reproduction instructions

Under the directory `post-augmentation-genotyping`, the command `snakemake --use-conda -j <number of cores>` can be used to run the entire pipeline.

Due to the large requirements of certain jobs, a high performance computing architecture will be needed. For running the pipeline in such an architecture, refer to the snakemake documentation (https://snakemake.readthedocs.io/en/v7.18.1/) for further details.
