# Code and Software Checklist

This README contains all the links of the repositories provided in the "Code Availability" section of the main text along with the README files' URL.

## Repository 1

The pipeline produces a merged, multi-sample VCF from the post-processed SVarp calls (Methods, Section “SVarp VCF postprocessing”).

URL of the Repository: https://github.com/eblerjana/long-read-1kg/tree/main/prepare-SVarp-callset

The URL of the README: https://github.com/eblerjana/long-read-1kg/blob/main/prepare-SVarp-callset/README.md


## Repository 2

The Graph Augmentation pipeline

URL of the Repository: https://github.com/eblerjana/long-read-1kg/tree/main/graph-augmentation-pipeline

The URL of the README: https://github.com/eblerjana/long-read-1kg/blob/main/graph-augmentation-pipeline/README.md


## Repository 3

The genotyping pipeline using Giggles (Methods, Section “Graph-aware genotyping with Giggles”), the associated preprocessing pipelines (Methods, Section “Preparing phased VCF panel”), SV calling with SVarp (Methods, Section “Variant calling using pangenome graph”), phasing pipeline (Methods, Section “Phasing with the ONT reads” and “Comparison of the WhatsHap phased VCFs against NYGC statistical phasing”), and the haplotype-tagging pipeline (Methods, Section “Haplotype-tagging of ONT reads”).

URL of the Repository: https://github.com/marschall-lab/project-ont-1kg/tree/main/snakemake_pipelines

The URL of the README:

 - Haplotype-tagging of ONT reads using WhatsHap for downstream tools: https://github.com/marschall-lab/project-ont-1kg/blob/main/snakemake_pipelines/haplotagging/README.md

 - Phasing experiments with WhatsHap for checking quality of reads: https://github.com/marschall-lab/project-ont-1kg/blob/main/snakemake_pipelines/phasing/README.md

 - Genotyping on the HPRC_mg_44+966 using Giggles along with preprocessing: https://github.com/marschall-lab/project-ont-1kg/blob/main/snakemake_pipelines/augmented_graph/README.md
 
 - Genotyping on the HPRC_mg using Giggles and calling variants with SVarp: https://github.com/marschall-lab/project-ont-1kg/blob/main/snakemake_pipelines/svarp-giggles_combined_pipeline/README.md


## Repository 4

Analysis script for genotyping.

URL of the Repository: https://github.com/marschall-lab/project-ont-1kg/tree/main/analysis_scripts/genotyping

The URL of the README: https://github.com/marschall-lab/project-ont-1kg/blob/main/analysis_scripts/README.md

## Repository 5

Further analysis scripts related to the base-calling, alignment to linear and graph-based references and SV calling.

URL of the Repository: https://github.com/1kg-ont-vienna/sv-analysis

The URL of the README: https://github.com/1kg-ont-vienna/sv-analysis/blob/main/README.md

## Repository 6

The structural variation annotation pipeline (Methods, Section “SVAN”).

URL of the Repository: https://github.com/REPBIO/SVAN

The URL of the README: ????

## Repository 7

The inversion analysis scripts for the simulations.

URL of the Repository: https://github.com/celiatsapalou/Simulations_ONT_Data

The URL of the README: https://github.com/celiatsapalou/Simulations_ONT_Data/blob/main/README.md

Pipeline for selecting regions for remapping and inversion calling.

URL of the Repository: https://github.com/celiatsapalou/Small_Inversions_Remap

The URL of the README: https://github.com/celiatsapalou/Small_Inversions_Remap/blob/master/README.md