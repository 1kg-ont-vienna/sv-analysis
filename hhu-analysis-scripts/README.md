# Details of the analysis scripts

The analysis scripts are for the genotypes that are the results of the `hhu-snakemake-pipelines/post-augmentation-genotyping` and `hhu-snakemake-pipelines/pre-augmentation-giggles-svarp` pipelines along with the phased genotypes created using ShapeIt5 by Dr Tobias Rausch.

Please refer to the `Makefile` for details on each section and the analysis done there.

## System Requirements

- python=3.9.18
- bcftools=1.19
- cyvcf2=0.30.22
- numpy=1.24.4
- scipy=1.8.1
- pandas=1.5.1
- matplotlib=3.5.2
- seaborn=0.12.2
- upsetplot=0.9.0
- bgzip=1.9.1
- awk=1.3.4

## Installation Guide

### Instructions

The *System requirements* given above can be installed using the command line using `apt` or `pip`. Please refer to the individual documentation for this.

## Instructions for use

### Input

Please refer to the `Makefile` about the inputs needed for running the scripts. Variables have been defined and comments have been provided of which files from the snakemake pipelines that need to be specified.

### Output

The main outputs of this pipeline are:

- Filtered VCFs using bcftools.
- QC tables for the VCFs containing various quality check parameters.
- Boxplots of SV lengths and SV counts
- Plots of log(allele count) vs log(number of variant sites)
- Quality assessments plots showing HWE agreements and comparison between allele frequencies in the panel and the genotypes.
- Audano SV growth curves
- Venn diagrams showing overlap between different consortiums
- Mendelian inconsistency calculation.

### Reproduction instructions

In the `Makefile`, there is a list of paths that need to be provided. These files come from the snakemake pipelines and after providing all the paths, the analysis can be run.

In the folder containing the `Makefile`, run `make all` for creating the output. To remove the outputs, run `make remove`.