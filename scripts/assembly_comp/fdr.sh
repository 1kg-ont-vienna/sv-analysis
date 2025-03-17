#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../mamba/bin:${PATH}

for INPUT in ${BASEDIR}/../release/giggles-genotyping/final-vcf.unphased.vcf.gz ${BASEDIR}/../release/giggles-genotyping/giggles-genotypes-biallelic.ac0-filtered.vcf.gz
do
    ID=`echo ${INPUT} | sed 's/^.*\///' | sed 's/.vcf.gz$//'`
    for SIZE in 50 250 SMALL
    do
	if [ ${SIZE} == "50" ]
	then
	    bcftools query -i '(STRLEN(ALT)>=(STRLEN(REF)+50))' -f "%ID\n" ${INPUT} > INS.ids
	    bcftools query -i '(STRLEN(REF)>=(STRLEN(ALT)+50))' -f "%ID\n" ${INPUT} > DEL.ids
	elif [ ${SIZE} == "250" ]
	then
	    bcftools query -i '(STRLEN(ALT)>=(STRLEN(REF)+250))' -f "%ID\n" ${INPUT} > INS.ids
	    bcftools query -i '(STRLEN(REF)>=(STRLEN(ALT)+250))' -f "%ID\n" ${INPUT} > DEL.ids
	else
	    bcftools query -i '(STRLEN(ALT)>=(STRLEN(REF)+50)) && (STRLEN(ALT)<(STRLEN(REF)+250))' -f "%ID\n" ${INPUT} > INS.ids
	    bcftools query -i '(STRLEN(REF)>=(STRLEN(ALT)+50)) && (STRLEN(REF)<(STRLEN(ALT)+250))' -f "%ID\n" ${INPUT} > DEL.ids
	fi
	for SV in INS DEL
	do
	    TP=`cat summary.${ID}_ALL.tsv | cut -f 1,6 | grep -w -Ff ${SV}.ids | awk '$2=="TP"' | wc -l`
	    FP=`cat summary.${ID}_ALL.tsv | cut -f 1,6 | grep -w -Ff ${SV}.ids | awk '$2=="FP"' | wc -l`
	    FDR=`echo "${FP} / ( ${TP} + ${FP} )" | bc -l`
	    echo ${ID} ${SV} ${SIZE} ${TP} ${FP} ${FDR}
	done
    done
done
