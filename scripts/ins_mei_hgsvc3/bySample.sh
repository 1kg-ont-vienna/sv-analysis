#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../mamba/bin:${PATH}

# Window for insertions
WIN=200
ONLY_CANONICAL=0

if [ ! -f MEI_Callset_T2T-CHM13.ALL.20240918.csv ]
then
    wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/release/Mobile_Elements/1.0/MEI_Callset_T2T-CHM13.ALL.20240918.csv.gz
    gunzip MEI_Callset_T2T-CHM13.ALL.20240918.csv.gz 
fi
if [ ! -f ont.vcf.gz ]
then
    bcftools view -O b -o tmp.bcf -S samples.tsv ${BASEDIR}/../release/giggles-genotyping/final-vcf.unphased.vcf.gz chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22
    bcftools index tmp.bcf
    bcftools annotate -c INFO -a ${BASEDIR}/../release/svan-annotation/final-vcf.unphased.SVAN_1.3.vcf.gz tmp.bcf | bgzip > ont.vcf.gz
    tabix ont.vcf.gz
    rm tmp.bcf tmp.bcf.csi
fi
echo "All MEIs"
bcftools query -i '(STRLEN(ALT)>=(STRLEN(REF)+50))' -f "%INFO/ITYPE_N\t%INFO/FAM_N\t%INFO/NOT_CANONICAL\n" ${BASEDIR}/../release/svan-annotation/final-vcf.unphased.SVAN_1.3.vcf.gz | egrep "partnered|solo" | egrep "Alu|L1|SVA" | cut -f 2 | sort | uniq -c
echo "Non-canonical MEIs"
bcftools query -i '(STRLEN(ALT)>=(STRLEN(REF)+50))' -f "%INFO/ITYPE_N\t%INFO/FAM_N\t%INFO/NOT_CANONICAL\n" ${BASEDIR}/../release/svan-annotation/final-vcf.unphased.SVAN_1.3.vcf.gz | egrep "partnered|solo" | egrep "Alu|L1|SVA" | cut -f 2- | sort | uniq -c

## HGSVC3
python3 convert.py -c MEI_Callset_T2T-CHM13.ALL.20240918.csv -w ${WIN} > hgsvc3.t2t.bed

## ONT
echo -e "mei\tsample\tdataset\tall\tshared\tunique\tfdr\tsensitivity" > stats.tsv
for SAMPLE in `cat samples.tsv`
do
    for TAG in Alu L1 SVA
    do
	if [ ${ONLY_CANONICAL} -eq 0 ]
	then
	    bcftools view -i '(STRLEN(ALT)>=(STRLEN(REF)+50))' -s ${SAMPLE} --min-ac 1 ont.vcf.gz | bcftools query -f "%CHROM\t%POS\t%ID\t%INFO/ITYPE_N\t%INFO/FAM_N\n" - | awk '{print $1"\t"($2-'${WIN}')"\t"($2+'${WIN}')"\t"$3"\t"$4"\t"$5;}' | egrep "partnered|solo" | grep -w "${TAG}" | sort -k1,1V -k2,2n | uniq > ont.${SAMPLE}.${TAG}.bed
	else
	    bcftools view -i '(STRLEN(ALT)>=(STRLEN(REF)+50)) && INFO/NOT_CANONICAL==0' -s ${SAMPLE} --min-ac 1 ont.vcf.gz | bcftools query -f "%CHROM\t%POS\t%ID\t%INFO/ITYPE_N\t%INFO/FAM_N\n" - | awk '{print $1"\t"($2-'${WIN}')"\t"($2+'${WIN}')"\t"$3"\t"$4"\t"$5;}' | egrep "partnered|solo" | grep -w "${TAG}" | sort -k1,1V -k2,2n | uniq > ont.${SAMPLE}.${TAG}.bed
	fi
	cat hgsvc3.t2t.bed | grep -w ${SAMPLE} | grep -w "${TAG}" | sort -k1,1V -k2,2n | uniq > hgsvc.${SAMPLE}.${TAG}.bed
	COUNTONT=`cat ont.${SAMPLE}.${TAG}.bed | cut -f 1-3 | sort | uniq | wc -l`
	COUNTHGSVC=`cat hgsvc.${SAMPLE}.${TAG}.bed | cut -f 1-3 | sort | uniq | wc -l`
	SHAREDONT=`bedtools intersect -a ont.${SAMPLE}.${TAG}.bed -b hgsvc.${SAMPLE}.${TAG}.bed -wao | awk '$7!="."' | cut -f 1-3 | sort | uniq | wc -l`
	UNIQONT=`bedtools intersect -a ont.${SAMPLE}.${TAG}.bed -b hgsvc.${SAMPLE}.${TAG}.bed -wao | awk '$7=="."' | cut -f 1-3 | sort | uniq | wc -l`
	SHAREDHGSVC=`bedtools intersect -a hgsvc.${SAMPLE}.${TAG}.bed -b ont.${SAMPLE}.${TAG}.bed -wao | awk '$6!="."' | cut -f 1-3 | sort | uniq | wc -l`
	UNIQHGSVC=`bedtools intersect -a hgsvc.${SAMPLE}.${TAG}.bed -b ont.${SAMPLE}.${TAG}.bed -wao | awk '$6=="."' | cut -f 1-3 | sort | uniq | wc -l`
	echo ${SAMPLE} ${TAG}
	FDRONT=`echo "${UNIQONT} / ${COUNTONT}" | bc -l`
	SENSONT=`echo "${SHAREDHGSVC} / ${COUNTHGSVC}" | bc -l`
	FDRHGSVC=`echo "${UNIQHGSVC} / ${COUNTHGSVC}" | bc -l`
	SENSHGSVC=`echo "${SHAREDONT} / ${COUNTONT}" | bc -l`
	echo -e "${TAG}\t${SAMPLE}\t1kG_ONT\t${COUNTONT}\t${SHAREDONT}\t${UNIQONT}\t${FDRONT}\t${SENSONT}" >> stats.tsv
	echo -e "${TAG}\t${SAMPLE}\tHGSVC3\t${COUNTHGSVC}\t${SHAREDHGSVC}\t${UNIQHGSVC}\t${FDRHGSVC}\t${SENSHGSVC}" >> stats.tsv
	rm ont.${SAMPLE}.${TAG}.bed hgsvc.${SAMPLE}.${TAG}.bed 
    done
done
rm hgsvc3.t2t.bed

Rscript bySample.R

for TAG in Alu L1 SVA
do
    # Avg. FDR
    echo "${TAG} average FDR:" `cat stats.tsv  | grep "ONT" | grep "${TAG}" | awk '{SUM+=$7;} END {print SUM/NR;}'`
    echo "${TAG} average sensitivity" `cat stats.tsv  | grep "ONT" | grep "${TAG}" | awk '{SUM+=$8;} END {print SUM/NR;}'`
done
