#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../mamba/bin:${PATH}

if [ ! -f 1KGP_3202.Illumina_ensemble_callset.freeze_V1.vcf.gz ]
then
    wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20210124.SV_Illumina_Integration/1KGP_3202.Illumina_ensemble_callset.freeze_V1.vcf.gz
    wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20210124.SV_Illumina_Integration/1KGP_3202.Illumina_ensemble_callset.freeze_V1.vcf.gz.tbi
fi

## Insertions
TOTAL=0
for TE in Alu L1 SVA
do
    COUNT=`bcftools query -i 'STRLEN(REF)<STRLEN(ALT)' -f "%ITYPE_N\t%FAM_N\n" ${BASEDIR}/../release/svan-annotation/final-vcf.unphased.SVAN_1.3.vcf.gz | awk '$2=="'${TE}'"' | wc -l`
    echo "Insertion" ${TE} ${COUNT}
    TOTAL=`expr ${TOTAL} + ${COUNT}`
done
echo "Non-reference MEI" ${TOTAL}
echo "Canonical/Non-canonical MEIs"
bcftools query -i 'STRLEN(REF)<STRLEN(ALT)' -f '%ITYPE_N\t%FAM_N\t%NOT_CANONICAL\n' ${BASEDIR}/../release/svan-annotation/final-vcf.unphased.SVAN_1.3.vcf.gz | grep -v "," | egrep -w "Alu|L1|SVA"  | cut -f 3 | sort | uniq -c

## Deletions
TOTAL=0
for TE in Alu L1 SVA
do
    COUNT=`bcftools query -i 'STRLEN(REF)>STRLEN(ALT)' -f "%DTYPE_N\t%FAM_N\n" ${BASEDIR}/../release/svan-annotation/final-vcf.unphased.SVAN_1.3.vcf.gz | awk '$2=="'${TE}'"' | wc -l`
    echo "Deletion" ${TE} ${COUNT}
    TOTAL=`expr ${TOTAL} + ${COUNT}`
done
echo "Reference MEI" ${TOTAL}

## Comparison to short-reads 
bcftools view ${BASEDIR}/../release/giggles-genotyping/final-vcf.unphased.vcf.gz | grep -m 1 "^#CHROM" | cut -f 10- | tr '\t' '\n' > samples.tsv
if [ ! -f ont.vcf.gz ]
then
    bcftools view -O b -o tmp.bcf --min-ac 1 -S samples.tsv ${BASEDIR}/../release/giggles-genotyping/final-vcf.unphased.vcf.gz chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX
    bcftools index tmp.bcf
    bcftools annotate -c INFO -a ${BASEDIR}/../release/svan-annotation/final-vcf.unphased.SVAN_1.3.vcf.gz tmp.bcf | bgzip > ont.vcf.gz
    tabix ont.vcf.gz
    rm tmp.bcf tmp.bcf.csi
fi
if [ ! -f nygc.vcf.gz ]
then
    bcftools view --min-ac 1 -S samples.tsv 1KGP_3202.Illumina_ensemble_callset.freeze_V1.vcf.gz chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX | bgzip > nygc.vcf.gz
    tabix nygc.vcf.gz
fi

## Insertions
for TE in Alu L1 SVA
do
    ONT=`bcftools query -i 'STRLEN(REF)<STRLEN(ALT)' -f "%ITYPE_N\t%FAM_N\n" ${BASEDIR}/../release/svan-annotation/final-vcf.unphased.SVAN_1.3.vcf.gz | awk '$2=="'${TE}'"' | wc -l`
    ALTTAG=`echo ${TE} | sed 's/Alu/<INS:ME:ALU>/' | sed 's/L1/<INS:ME:LINE1>/' | sed 's/SVA/<INS:ME:SVA>/'`
    NYGC=`zcat nygc.vcf.gz | grep -v "^#" | cut -f 5 | grep -w "${ALTTAG}" | wc -l`
    INC=`echo " ( ( ${ONT} - ${NYGC} )  / ${NYGC} ) "| bc -l`
    echo "Comparison to short-reads:" ${TE} ${INC}
done
