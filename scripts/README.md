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

The scripts assume these release data sets are available in the subfolder `./release/`.

## Alignment statistics

To calculate the alignment error rate, genome coverage and other statistics.

`cd qc/ && ./qc.sh`

## MEI statistics

The MEI statistics can be computed using

`cd mei_statistics/ && ./stats.sh`

## MEI allele length distribution

The comparison of sequence-resolved transposable element insertions from ONT compared to illumina-based MEI predictions

`cd mei_allele_length/ && ./len.sh`

## VNTR density by Genome in a Bottle genomic stratifications

The analysis of VNTR density can be reproduced using

`cd vntr_density/ && ./vntr.sh`

## Comparison of the SAGA SV call set to HGSVC3

To run the call set comparison for the 16 samples shared between the 2 projects use:

`cd assembly_comp/ && ./hgsvc3.sh`

Further scripts are included to reproduce the false negative density plot (`./fn.plot.sh`), compute the FDRs (`./fdr.sh`) and run the sample-specific FDR analysis (`./fdrBySample.sh`). The assembly comparison revealed 556 SVs that are present in every sample of the 1000 Genomes ONT study that are likely representing misalignments or graph construction errors (`graph.errors.svs.gz`). The TRF annotation track for CHM13 (T2T) can be downloaded from the UCSC table browser (`repeat.trf.chm13.msk.gz`).

## Comparing SAGA MEIs to multi-platform whole-genome assemblies from HGSVC3

The MEI comparison to HGSVC3 can be recomputed using:

`cd ins_mei_hgsvc3/ && ./bySample.sh`

For transposable element deletions from the CHM13 reference:

`cd del_mei_hgsvc3/ && ./bySample.sh`

## Intensity rank sum (IRS) test for deletions

The intensity rank sum test using Affymetrix SNP6 arrays.

`cd irs/ && ./runIRS.sh`

## Lifting of SVs to GRCh38

This workflow uses the phased SVs for each sample to implant SVs into CHM13 haplotypes, which are then aligned to GRCh38 and called via [svim-asm](https://github.com/eldariont/svim-asm).

`cd svim-asm/ && ./svim.sh NA20504`

All by-sample VCFs are then merged into a multi-sample VCF using:

```
bcftools merge --no-version -0 -m snp-ins-del -O b -o merged.bcf *.norm.bcf
bcftools index merged.bcf
```

To remove SVs that are present in all 908 samples, i.e. CHM13 differences compared to GRCh38, we use bcftools: 

```
bcftools view merged.bcf chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX | bcftools +fill-tags - -- -t all | awk 'BEGIN {FS="\t"; IFS="\t"; OFS="\t";} { if (substr($1,1,1)!="#") { $3=sprintf("SvimAsm%08d", NR); }; print $0;}' | grep -v "AC=1816" | bcftools view -O b -o svim.asm.hg38.bcf --min-ac 1 -i '(STRLEN(REF)>=(STRLEN(ALT)+50)) || (STRLEN(ALT)>=(STRLEN(REF)+50))' -`
bcftools index svim.asm.hg38.bcf
```

## Comparison of SVs to gnomAD

Using the lifted GRCh38 call set, a comparison to gnomAD can be carried out using:

`cd gnomad/ && ./novelty.sh`


## Credits

[HTSlib](https://github.com/samtools/htslib), [SAMtools](https://github.com/samtools/samtools), [BCFtools](https://github.com/samtools/bcftools) and [BEDtools](https://github.com/arq5x/bedtools2) for alignment, VCF and BED file processing. [Minimap2](https://github.com/lh3/minimap2) for long-read alignment to linear reference genomes, [Minigraph](https://github.com/lh3/minigraph) for pangenome graph construction and alignment. [ggplot2](https://ggplot2.tidyverse.org/index.html) for plotting.

