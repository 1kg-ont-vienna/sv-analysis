# Code and Software Checklist

This README contains all the links of the repositories provided in the "Code Availability" section of the main text along with the README files' URL.


## Repository 1

Repository: https://github.com/eblerjana/long-read-1kg

Commit: 44a1752

Description: Pipeline to create a multi-sample SVarp VCF from single-sample SVarp calls (used as one of the callsets going into graph augmentation) and the graph augmentation pipeline.


## Repository 2

Repository: https://github.com/samarendra-pani/giggles

Commit: 5226884 (version 1.0)

Description: The SV genotyper Giggles which uses long read and pangenome reference.


## Repository 3

Repository: https://github.com/marschall-lab/project-ont-1kg

Commit: 32d896d

Description:

The following snakemake pipelines are hosted in this repository:

- Haplotagging of Aligned Reads under `haplotagging`: The pipeline tags the aligned reads of this study as originating from haplotype 1 or 2 using `whatshap`.
- Phasing Experiments under `phasing`: The pipeline phases the NYGC raw genotyes using the aligned reads with `whatshap` and we perform QC using the NYGC statistical phased VCF.
- Running Giggles and SVarp on the HPRC_mg graph under `pre-augmentation-giggles-svarp`: The pipeline pre-processes and does SV discovery with `SVarp` and SV genotyping with `giggles` using HPRC_mg.
- Running Giggles on the HPRC_mg_44+966 graph under `post-augmentation-genotyping`: The pipeline pre-processes and does SV genotyping with `giggles` using HPRC_mg_44+966.
- Running QCs based on coverage and read N50 stratification under `genotype-sample-subet-analysis`: The pipeline bins the samples based on coverage and read N50 and runs QC to check the effects of these two variables on genotyping and SV discovery.
- Investigating recovery of added SVs during genotyping under `sv-recovery`: The pipeline investigates the SVs that were introduced into the graph to create HPRC_mg_44+966 and whether they were retrieved during genotyping.
- Annotating the alleles in HPRC_mg_44+966 based on the ancestral allele found in Chimpanzee under `chimpanzee-analysis`: The pipeline identifies ancestral allele in the chimpanzee using alignment of the chimpanzee assembly to HPRC_mg_44+966 and tags the VCFs produced from `giggles`.
- VNTR genotyping of the cohort using vamos under `vamos-analysis`: The pipeline runs the VNTR calling using `vamos` and subsequent QC with the HGSVC3 assemblies.

The repository also hosts analysis scripts for the genotypes and SVAN annotations under `analysis_scripts`. The repository also has the curation work for VNTR-associated diseases under `vntr-based-rare-disease-curation`.


## Repository 4

Repository: https://github.com/REPBIO/SVAN

Commit: XXX

Description: The SV Annotation (SVAN) scripts.


## Repository 5

Repository: https://github.com/carstenhain/SV_homology

Commit: 534b216

Description: Scripts for finding and annotating homologus flanks at SV breakpoints.


## Repository 6

Repository: https://github.com/celiatsapalou/Simulations_ONT_Data

Commit: cf0513b

Description: Pipeline to simulate ONT data on a genome with small inversions for benchmarking genome mappers and callers.


## Repository 7

Repository: https://github.com/celiatsapalou/Small_Inversions_Remap

Commit: f9d8a44

Description: Pipeline to export high-mismatch regions, remap with NGMLR, and call inversions.


## Repository 8

Repository: https://github.com/RMoreiraP/GeONTIpe

Commit: 1b5db07

Description: Pipeline to genotype inversions using ONT reads.


## Repository 9

Respository: https://github.com/hugocarmaga/rare-disease-analysis

Commit: add69e9

Pipeline to filter and get some statistics on rare disease samples' SVs after comparing them with big SV resources.


## Repository 10

Repository: https://github.com/hugocarmaga/variant-calling

Commit: XXX

Description: SV calling pipeline from HiFi reads to structural variation VCF files