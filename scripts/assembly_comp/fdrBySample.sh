#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../mamba/bin:${PATH}

echo -e "sample\tvcf\tsvtype\tclass\ttp\tfp\tfdr\ttrf_svs\ttrf_fraction" > fdr.by_sample.tsv
for INPUT in ${BASEDIR}/../release/giggles-genotyping/final-vcf.unphased.vcf.gz
do
    ID=`echo ${INPUT} | sed 's/^.*\///' | sed 's/.vcf.gz$//'`

    ## Overall TRF counts
    bcftools query -i '(STRLEN(ALT)>=(STRLEN(REF)+250)) || (STRLEN(REF)>=(STRLEN(ALT)+250))' -f "%CHROM\t%POS\t%ID\n" ${INPUT} | awk '{print $1"\t"$2"\t"($2+1)"\t"$3;}' > svs.bed
    TOTAL=`cat svs.bed | wc -l`
    TRFSV=`bedtools intersect -a svs.bed -b <(zcat repeat.trf.chm13.msk.gz | tail -n +2 | awk '$5*$6>=50')  | sort | uniq | wc -l`
    TRFFRAC=`echo "${TRFSV} / ${TOTAL}" | bc -l`
    echo "TRF fraction for SVs >=250bp" ${TRFSV} ${TOTAL} ${TRFFRAC}
    bcftools query -i '((STRLEN(ALT)>=(STRLEN(REF)+50)) && (STRLEN(ALT)<(STRLEN(REF)+250))) || ((STRLEN(REF)>=(STRLEN(ALT)+50)) && (STRLEN(REF)<(STRLEN(ALT)+250)))' -f "%CHROM\t%POS\t%ID\n" ${INPUT} | awk '{print $1"\t"$2"\t"($2+1)"\t"$3;}' > svs.bed
    TOTAL=`cat svs.bed | wc -l`
    TRFSV=`bedtools intersect -a svs.bed -b <(zcat repeat.trf.chm13.msk.gz | tail -n +2 | awk '$5*$6>=50')  | sort | uniq | wc -l`
    TRFFRAC=`echo "${TRFSV} / ${TOTAL}" | bc -l`
    echo "TRF fraction for SVs <250bp" ${TRFSV} ${TOTAL} ${TRFFRAC}
    rm svs.bed
    
    ## TRF overlap
    bcftools query -f "%CHROM\t%POS\t%ID\n" ${INPUT} | awk '{print $1"\t"$2"\t"($2+1)"\t"$3;}' > svs.bed
    bedtools intersect -a svs.bed -b <(zcat repeat.trf.chm13.msk.gz | tail -n +2 | awk '$5*$6>=50')   | cut -f 4 | sort | uniq > trf.svs
    for SUMTSV in summary.${ID}*_[HN][GA][0-9]*.tsv
    do
	SAMPLE=`echo ${SUMTSV} | sed 's/.tsv$//' | sed 's/^.*_//'`
	echo ${SAMPLE}
	for SIZE in SMALL LARGE
	do
	    if [ ${SIZE} == "LARGE" ]
	    then
		bcftools query -i '(STRLEN(ALT)>=(STRLEN(REF)+250))' -f "%ID\n" ${INPUT} > INS.ids
		bcftools query -i '(STRLEN(REF)>=(STRLEN(ALT)+250))' -f "%ID\n" ${INPUT} > DEL.ids
	    else
		bcftools query -i '(STRLEN(ALT)>=(STRLEN(REF)+50)) && (STRLEN(ALT)<(STRLEN(REF)+250))' -f "%ID\n" ${INPUT} > INS.ids
		bcftools query -i '(STRLEN(REF)>=(STRLEN(ALT)+50)) && (STRLEN(REF)<(STRLEN(ALT)+250))' -f "%ID\n" ${INPUT} > DEL.ids
	    fi
	    for SV in INS DEL
	    do
		TP=`cat ${SUMTSV} | cut -f 1,6 | grep -w -Ff ${SV}.ids | awk '$2=="TP"' | wc -l`
		FP=`cat ${SUMTSV} | cut -f 1,6 | grep -w -Ff ${SV}.ids | awk '$2=="FP"' | wc -l`
		FDR=`echo "${FP} / ( ${TP} + ${FP} )" | bc -l`
		TRFSV=`cat ${SUMTSV} | cut -f 1,6 | grep -w -Ff ${SV}.ids | grep -w -Ff trf.svs | awk '$2=="TP" || $2=="FP"' | wc -l`
		FRACTRF=`echo "${TRFSV} / ( ${TP} + ${FP} )" | bc -l`
		echo -e "${SAMPLE}\t${ID}\t${SV}\t${SIZE}\t${TP}\t${FP}\t${FDR}\t${TRFSV}\t${FRACTRF}" >> fdr.by_sample.tsv
	    done
	    rm DEL.ids INS.ids
	done
    done
    rm trf.svs svs.bed
done
Rscript fdrBySample.R
