#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../mamba/bin:${PATH}

if [ ! -f 1KGP_3202.Illumina_ensemble_callset.freeze_V1.vcf.gz ]
then
    wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20210124.SV_Illumina_Integration/1KGP_3202.Illumina_ensemble_callset.freeze_V1.vcf.gz
fi

if [ ! -f nygc.mei ]
then
    bcftools view 1KGP_3202.Illumina_ensemble_callset.freeze_V1.vcf.gz chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22  | bcftools query -f "%CHROM\t%POS\t%SVLEN\t%ALT\n" - | grep -v "<INS:ME>" |  grep "INS:ME" | sed 's/<INS:ME://' | sed 's/>//' > nygc.mei
    cat nygc.mei  | sed 's/ALU/Alu/' | sed 's/LINE1/L1/' > nygc.mei.tmp && mv nygc.mei.tmp nygc.mei
fi

if [ ! -f ont.vcf.gz ]
then
    bcftools view -O b -o tmp.bcf ${BASEDIR}/../release/giggles-genotyping/final-vcf.unphased.vcf.gz chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22
    bcftools index tmp.bcf
    bcftools annotate -c INFO -a ${BASEDIR}/../release/svan-annotation/final-vcf.unphased.SVAN_1.3.vcf.gz tmp.bcf | bgzip > ont.vcf.gz
    tabix ont.vcf.gz
    rm tmp.bcf tmp.bcf.csi
    bcftools query -f "%CHROM\t%POS\t%INS_LEN\t%ITYPE_N\t%FAM_N\n" ont.vcf.gz | egrep "Alu|L1|SVA" | cut -f 1,2,3,5 | grep -v "," | awk '$3!="."'  > ont.mei
fi
