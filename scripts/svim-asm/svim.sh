#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../mamba/bin:${PATH}

REF=${BASEDIR}/../genome/t2t.fa
HG38=${BASEDIR}/../genome/hg38.fa
export REF_CACHE=${BASEDIR}/../genome/cache_t2t/%2s/%2s/%s

THREADS=8

SAMPLE=${1}

echo ${SAMPLE}
mkdir -p ${SAMPLE}
cd ${SAMPLE}
# Create assemblies
bcftools view -O b -o ${SAMPLE}.bcf -s ${SAMPLE} --min-ac 1 ${BASEDIR}/../release/shapeit4-phased-callset/final-vcf.phased.vcf.gz
bcftools index ${SAMPLE}.bcf
cat ${REF} | bcftools consensus -s ${SAMPLE} -H 1 ${SAMPLE}.bcf > ${SAMPLE}.h1.fa
cat ${REF} | bcftools consensus -s ${SAMPLE} -H 2 ${SAMPLE}.bcf > ${SAMPLE}.h2.fa

# Map
minimap2 -a -x asm5 --cs -r2k -t 24 ${HG38} ${SAMPLE}.h1.fa | samtools sort -m 4G -@ 8 -o ${SAMPLE}.h1.bam -
samtools index ${SAMPLE}.h1.bam
minimap2 -a -x asm5 --cs -r2k -t 24 ${HG38} ${SAMPLE}.h2.fa | samtools sort -m 4G -@ 8 -o ${SAMPLE}.h2.bam -
samtools index ${SAMPLE}.h2.bam

# Call SVs
svim-asm diploid ${SAMPLE} ${SAMPLE}.h1.bam ${SAMPLE}.h2.bam ${HG38}
cd ../
bcftools view -O b -o ${SAMPLE}.bcf ${SAMPLE}/${SAMPLE}/variants.vcf
bcftools index ${SAMPLE}.bcf
rm -rf ${SAMPLE}/

## Normalize
bcftools view --min-ac 1 -i '(STRLEN(REF)>=(STRLEN(ALT)+50)) || (STRLEN(ALT)>=(STRLEN(REF)+50))' ${SAMPLE}.bcf | sed "s/\tSample$/\t${SAMPLE}/" | bcftools norm -c w -a -f ${HG38} -m -any --no-version - | bcftools view -O b -o ${SAMPLE}.norm.bcf --min-ac 1 -i '(STRLEN(REF)>=(STRLEN(ALT)+50)) || (STRLEN(ALT)>=(STRLEN(REF)+50))' - 
bcftools index ${SAMPLE}.norm.bcf
