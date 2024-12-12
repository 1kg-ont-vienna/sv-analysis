#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../mamba/bin:${PATH}

## Insertions
TOTAL=0
for TE in Alu L1 SVA
do
    COUNT=`bcftools query -i 'STRLEN(REF)<STRLEN(ALT)' -f "%ITYPE_N\t%FAM_N\n" ../release/svan-annotation/final-vcf.unphased.SVAN_1.3.vcf.gz | awk '$2=="'${TE}'"' | wc -l`
    echo "Insertion" ${TE} ${COUNT}
    TOTAL=`expr ${TOTAL} + ${COUNT}`
done
echo "Non-reference MEI" ${TOTAL}

## Deletions
TOTAL=0
for TE in Alu L1 SVA
do
    COUNT=`bcftools query -i 'STRLEN(REF)>STRLEN(ALT)' -f "%DTYPE_N\t%FAM_N\n" ../release/svan-annotation/final-vcf.unphased.SVAN_1.3.vcf.gz | awk '$2=="'${TE}'"' | wc -l`
    echo "Deletion" ${TE} ${COUNT}
    TOTAL=`expr ${TOTAL} + ${COUNT}`
done
echo "Reference MEI" ${TOTAL}

