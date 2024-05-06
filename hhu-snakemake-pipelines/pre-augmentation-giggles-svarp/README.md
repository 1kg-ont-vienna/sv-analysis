# Details of the Pipeline

- This snakemake pipeline is to genotype and call variants using the longread dataset and the HPRC_mg graph.
- Genotypes are inferred with Giggles and SVs are called with SVarp.
- The SV calls of SVarp from this pipeline are then input to the graph augmentation pipeline to create the HPRC_mg_44+966 graph.
- The Giggles genotypes here act as a benchmark for comparison to the genotypes on HPRC_mg_44+966.

## System Requirements

- python=3.7.12
- pandas=1.3.5
- snakemake=7.18.1
- conda=22.9.0
- bgzip=1.9
- tabix=1.9
- minigraph=0.20-r559
- gaftools (commit ID: [903921ccc7ef4ffb6dd150974c080e75b5d72aac](https://github.com/marschall-lab/gaftools/tree/903921ccc7ef4ffb6dd150974c080e75b5d72aac))
- wtdbg2=2.5
- svarp (commit ID: [0acba75ebdfdd292d57e1bd133d852f6371ab677](https://github.com/asylvz/SVarp/tree/0acba75ebdfdd292d57e1bd133d852f6371ab677))
- singularity=3.5.2
- giggles=1.0
- pav=2.2.4.1

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

- svarp
    - link to github: https://github.com/asylvz/SVarp
    - the path to the executable needs to be provided after installation

- giggles
    - link to github: https://github.com/samarendra-pani/giggles
    - for this pipeline, giggles needs to be installed in a conda environment called `giggles-env`. Instructions on conda installation is provided in the github wiki.

- wtdbg2
    - link to github: https://github.com/ruanjue/wtdbg2
    - the repository has to be locally cloned, installed and the path to the directory needs to be provided.

- pav
    - https://github.com/EichlerLab/pav
    - the pipeline uses a singularity container of PAV.
    - The specific version used for the pipeline is not available online but can be obtained from Samarendra Pani. Contact at samarendra.pani@hhu.de

## Instructions for use

### Input

Running the pipeline on our data will require downloading the entire dataset. Refer to the *config.yaml* file where it states the required data for the pipeline and where to download them. The required directories and files need to be specified in the variables and then the pipeline can be run.

### Output

The main outputs of this pipeline are:

- sorted GAF reads at `result/svarp-giggles/data/gaf/<sample>.sorted.gaf.gz`
- the VCF panel used as input for giggles at `result/svarp-giggles/data/vcf/panel-multiallelic.vcf`
- giggles genotypes in multiallelic (`result/svarp-giggles/genotypes/multisample-multiallelic.vcf.gz`) and biallelic form (`result/svarp-giggles/genotypes/multisample-biallelic.vcf.gz`)
- SVarp svtigs at `result/svarp-giggles/svarp/<sample>/`
- SVarp svtigs processed using svim-asm and pav at `result/svarp-giggles/svarp/<sample>/pav_<ref>/pav_svtigs_merged.vcf` for both T2T and hg38 reference.

### Reproduction instructions

Under the directory `pre-augmentation-giggles-svarp`, the command `snakemake --use-conda -j <number of cores> all` can be used to run the entire pipeline.

Due to the large requirements of certain jobs, a high performance computing architecture will be needed. For running the pipeline in such an architecture, refer to the snakemake documentation (https://snakemake.readthedocs.io/en/v7.18.1/) for further details.
