#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../mamba/bin:${PATH}

# HGSVC3
if [ ! -f hgsvc3.vcf.gz ]
then
    wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/release/Variant_Calls/1.0/T2T-CHM13/variants_T2T-CHM13_sv_insdel_alt_HGSVC2024v1.0.vcf.gz
    tabix variants_T2T-CHM13_sv_insdel_alt_HGSVC2024v1.0.vcf.gz
    bcftools view variants_T2T-CHM13_sv_insdel_alt_HGSVC2024v1.0.vcf.gz chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 | grep -v '##contig=<ID=chrM' | bcftools +fill-tags - -- -t all | bgzip > hgsvc3.vcf.gz
    tabix hgsvc3.vcf.gz
fi

# ONT
for INPUT in ${BASEDIR}/../release/giggles-genotyping/final-vcf.unphased.vcf.gz ${BASEDIR}/../release/giggles-genotyping/giggles-genotypes-biallelic.ac0-filtered.vcf.gz
do
    ID=`echo ${INPUT} | sed 's/^.*\///' | sed 's/.vcf.gz$//'`
    echo ${ID}
    
    ## Fetch autosomes
    rm -rf ont.vcf.gz*
    bcftools view ${INPUT} chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 | bcftools +fill-tags - -- -t all | bgzip > ont.vcf.gz
    tabix ont.vcf.gz
    
    ## Collect SV properties
    bcftools query -i '(STRLEN(ALT)>=(STRLEN(REF)+50)) || (STRLEN(REF)>=(STRLEN(ALT)+50))' -f "%ID\t%MAF\t%AC" hgsvc3.vcf.gz | sed 's/-/\t/g' | awk '{print $1"-"$2"-"$3"-"$4"\t"$3"\t"$4"\t"$5"\t"$6;}' > sv.${ID}.properties
    bcftools query -i '(STRLEN(ALT)>=(STRLEN(REF)+50))' -f "%ID\t%MAF\t%AC" ont.vcf.gz | sed 's/-/\t/g' | awk '{print $1"-"$2"-"$3"-"$4"-"$5"\tINS\t"$5"\t"$6"\t"$7;}' >> sv.${ID}.properties
    bcftools query -i '(STRLEN(REF)>=(STRLEN(ALT)+50))' -f "%ID\t%MAF\t%AC" ont.vcf.gz | sed 's/-/\t/g' | awk '{print $1"-"$2"-"$3"-"$4"-"$5"\tDEL\t"$5"\t"$6"\t"$7;}' >> sv.${ID}.properties
    
    # Intersect to common samples
    bcftools view hgsvc3.vcf.gz | grep -m 1 "^#CHROM" | cut -f 10- | tr '\t' '\n' | sort | uniq > hgsvc3.samples
    bcftools view ont.vcf.gz | grep -m 1 "^#CHROM" | cut -f 10- | tr '\t' '\n' | sort | uniq > ont.samples
    sort hgsvc3.samples ont.samples | sort | uniq -d > common.samples
    rm hgsvc3.samples ont.samples
    rm -f hgsvc3.common.vcf.gz*
    
    # Iterate samples
    for SAMPLE in ALL `cat common.samples`
    do
	echo ${SAMPLE}
	if [ ${SAMPLE} == "ALL" ]
	then
	    bcftools view --min-ac 1 -S common.samples hgsvc3.vcf.gz | bgzip > hgsvc3.common.vcf.gz
	else
	    bcftools view --min-ac 1 -s ${SAMPLE} hgsvc3.vcf.gz | bgzip > hgsvc3.common.vcf.gz
	fi
	tabix hgsvc3.common.vcf.gz
	rm -f ont.common.vcf.gz*
	if [ ${SAMPLE} == "ALL" ]
	then
	    bcftools view --min-ac 1 -S common.samples ont.vcf.gz | bgzip > ont.common.vcf.gz
	else
	    bcftools view --min-ac 1 -s ${SAMPLE} ont.vcf.gz | bgzip > ont.common.vcf.gz
	fi
	tabix ont.common.vcf.gz

	## Truvari
	rm -rf ${ID}_${SAMPLE}/
	truvari bench --pick multi -p 0.5 -P 0.25 -b hgsvc3.common.vcf.gz -c ont.common.vcf.gz -f ${BASEDIR}/../genome/t2t.fa -o ${ID}_${SAMPLE}/

	## Collect summary stats
	cat sv.${ID}.properties | grep -w -Ff <(bcftools query -i '(STRLEN(ALT)>=(STRLEN(REF)+50)) || (STRLEN(REF)>=(STRLEN(ALT)+50))' -f "%ID\n" ${ID}_${SAMPLE}/tp-comp.vcf.gz) | sed 's/$/\tTP/' > summary.${ID}_${SAMPLE}.tsv
	cat sv.${ID}.properties	| grep -w -Ff <(bcftools query -i '(STRLEN(ALT)>=(STRLEN(REF)+50)) || (STRLEN(REF)>=(STRLEN(ALT)+50))' -f "%ID\n" ${ID}_${SAMPLE}/fp.vcf.gz) | sed 's/$/\tFP/' >> summary.${ID}_${SAMPLE}.tsv
	cat sv.${ID}.properties | grep -w -Ff <(bcftools query -i '(STRLEN(ALT)>=(STRLEN(REF)+50)) || (STRLEN(REF)>=(STRLEN(ALT)+50))' -f "%ID\n" ${ID}_${SAMPLE}/fn.vcf.gz) | sed 's/$/\tFN/' >> summary.${ID}_${SAMPLE}.tsv
	Rscript hgsvc3.R summary.${ID}_${SAMPLE}.tsv
	
	## Clean-up
	rm -rf refine.base.vcf.gz* refine.comp.vcf.gz* ont.common.vcf.gz* hgsvc3.common.vcf.gz*
    done
    rm -rf ont.vcf.gz*
done
