#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../mamba/bin:${PATH}

for INPUT in ${BASEDIR}/../release/giggles-genotyping/final-vcf.unphased.vcf.gz ${BASEDIR}/../release/giggles-genotyping/giggles-genotypes-biallelic.ac0-filtered.vcf.gz
do
    ID=`echo ${INPUT} | sed 's/^.*\///' | sed 's/.vcf.gz$//'`
    for SIZE in all small large
    do
	if [ ${SIZE} == "all" ]
	then
	    bcftools query -i '(STRLEN(ALT)>=(STRLEN(REF)+50))' -f "%ID\n" ${INPUT} > INS.ids
	    bcftools query -i '(STRLEN(REF)>=(STRLEN(ALT)+50))' -f "%ID\n" ${INPUT} > DEL.ids
	    bcftools query -i '(STRLEN(ALT)>=(STRLEN(REF)+50))' -f "%ID\n" ${ID}_ALL/fn.vcf.gz > INS_FN.ids
	    bcftools query -i '(STRLEN(REF)>=(STRLEN(ALT)+50))' -f "%ID\n" ${ID}_ALL/fn.vcf.gz > DEL_FN.ids
	elif [ ${SIZE} == "large" ]
	then
	    bcftools query -i '(STRLEN(ALT)>=(STRLEN(REF)+250))' -f "%ID\n" ${INPUT} > INS.ids
	    bcftools query -i '(STRLEN(REF)>=(STRLEN(ALT)+250))' -f "%ID\n" ${INPUT} > DEL.ids
	    bcftools query -i '(STRLEN(ALT)>=(STRLEN(REF)+250))' -f "%ID\n" ${ID}_ALL/fn.vcf.gz > INS_FN.ids
	    bcftools query -i '(STRLEN(REF)>=(STRLEN(ALT)+250))' -f "%ID\n" ${ID}_ALL/fn.vcf.gz > DEL_FN.ids
	else
	    bcftools query -i '(STRLEN(ALT)>=(STRLEN(REF)+50)) && (STRLEN(ALT)<(STRLEN(REF)+250))' -f "%ID\n" ${INPUT} > INS.ids
	    bcftools query -i '(STRLEN(REF)>=(STRLEN(ALT)+50)) && (STRLEN(REF)<(STRLEN(ALT)+250))' -f "%ID\n" ${INPUT} > DEL.ids
	    bcftools query -i '(STRLEN(ALT)>=(STRLEN(REF)+50)) && (STRLEN(ALT)<(STRLEN(REF)+250))' -f "%ID\n" ${ID}_ALL/fn.vcf.gz > INS_FN.ids
	    bcftools query -i '(STRLEN(REF)>=(STRLEN(ALT)+50)) && (STRLEN(REF)<(STRLEN(ALT)+250))' -f "%ID\n" ${ID}_ALL/fn.vcf.gz > DEL_FN.ids
	fi
	for SV in INS DEL
	do
	    TP=`cat summary.${ID}_ALL.tsv | cut -f 1,6 | grep -w -Ff ${SV}.ids | awk '$2=="TP"' | wc -l`
	    FP=`cat summary.${ID}_ALL.tsv | cut -f 1,6 | grep -w -Ff ${SV}.ids | awk '$2=="FP"' | wc -l`
	    FN=`cat summary.${ID}_ALL.tsv | cut -f 1,6 | grep -w -Ff ${SV}_FN.ids | awk '$2=="FN"' | wc -l`
	    FDR=`echo "${FP} / ( ${TP} + ${FP} )" | bc -l`
	    TPR=`echo "${TP} / ( ${TP} + ${FN} )" | bc -l`
	    echo ${ID} ${SV} ${SIZE} ${TP} ${FP} ${FN} ${FDR} ${TPR}
	done
	rm INS.ids DEL.ids INS_FN.ids DEL_FN.ids
    done    
    for MAF in all rare common
    do
	bcftools query -i '(STRLEN(ALT)>=(STRLEN(REF)+50))' -f "%ID\n" ${INPUT} > INS.ids
	bcftools query -i '(STRLEN(REF)>=(STRLEN(ALT)+50))' -f "%ID\n" ${INPUT} > DEL.ids
	bcftools query -i '(STRLEN(ALT)>=(STRLEN(REF)+50))' -f "%ID\n" ${ID}_ALL/fn.vcf.gz > INS_FN.ids
	bcftools query -i '(STRLEN(REF)>=(STRLEN(ALT)+50))' -f "%ID\n" ${ID}_ALL/fn.vcf.gz > DEL_FN.ids
	if [ ${MAF} == "rare" ]
	then
	    for F in *.ids
	    do
		cat ${F} | grep -w -Ff <(cat sv.${ID}.properties  | awk '$4<0.1' | cut -f 1 | sort | uniq) > ${F}.tmp
		mv ${F}.tmp ${F}
	    done
	elif [ ${MAF} == "common" ]
	then
	    for F in *.ids
	    do
		cat ${F} | grep -w -Ff <(cat sv.${ID}.properties  | awk '$4>=0.1' | cut -f 1 | sort | uniq) > ${F}.tmp
		mv ${F}.tmp ${F}
	    done
	fi
	for SV in INS DEL
	do
	    TP=`cat summary.${ID}_ALL.tsv | cut -f 1,6 | grep -w -Ff ${SV}.ids | awk '$2=="TP"' | wc -l`
	    FP=`cat summary.${ID}_ALL.tsv | cut -f 1,6 | grep -w -Ff ${SV}.ids | awk '$2=="FP"' | wc -l`
	    FN=`cat summary.${ID}_ALL.tsv | cut -f 1,6 | grep -w -Ff ${SV}_FN.ids | awk '$2=="FN"' | wc -l`
	    FDR=`echo "${FP} / ( ${TP} + ${FP} )" | bc -l`
	    TPR=`echo "${TP} / ( ${TP} + ${FN} )" | bc -l`
	    echo ${ID} ${SV} ${MAF} ${TP} ${FP} ${FN} ${FDR} ${TPR}
	done
	rm INS.ids DEL.ids INS_FN.ids DEL_FN.ids
    done
done
