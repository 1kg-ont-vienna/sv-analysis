#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../mamba/bin:${PATH}

if [ ! -f 1KGP_3202.Illumina_ensemble_callset.freeze_V1.vcf.gz ]
then
    wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20210124.SV_Illumina_Integration/1KGP_3202.Illumina_ensemble_callset.freeze_V1.vcf.gz
    wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20210124.SV_Illumina_Integration/1KGP_3202.Illumina_ensemble_callset.freeze_V1.vcf.gz.tbi
fi

## Comparison to short-reads 
bcftools view ${BASEDIR}/../release/giggles-genotyping/final-vcf.unphased.vcf.gz | grep -m 1 "^#CHROM" | cut -f 10- | tr '\t' '\n' > samples.tsv
exit;


## Insertions
TOTAL=0
for TE in Alu L1 SVA
do
    COUNT=`bcftools query -i 'STRLEN(REF)<STRLEN(ALT)' -f "%ITYPE_N\t%FAM_N\n" ${BASEDIR}/../release/svan-annotation/final-vcf.unphased.SVAN_1.3.vcf.gz | awk '$2=="'${TE}'"' | wc -l`
    echo "Insertion" ${TE} ${COUNT}
    TOTAL=`expr ${TOTAL} + ${COUNT}`
done
echo "Non-reference MEI" ${TOTAL}

## Deletions
TOTAL=0
for TE in Alu L1 SVA
do
    COUNT=`bcftools query -i 'STRLEN(REF)>STRLEN(ALT)' -f "%DTYPE_N\t%FAM_N\n" ${BASEDIR}/../release/svan-annotation/final-vcf.unphased.SVAN_1.3.vcf.gz | awk '$2=="'${TE}'"' | wc -l`
    echo "Deletion" ${TE} ${COUNT}
    TOTAL=`expr ${TOTAL} + ${COUNT}`
done
echo "Reference MEI" ${TOTAL}

