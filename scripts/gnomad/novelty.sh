#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../mamba/bin:${PATH}

if [ ! -f gnomad.v4.1.sv.sites.bed.gz ]
then
    wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.sites.bed.gz
fi
if [ ! -f gnomad.ins.bed ]
then
    zcat gnomad.v4.1.sv.sites.bed.gz | cut -f 1-5 | grep "INS" | sort -k1,1V -k2,2n | uniq > gnomad.ins.bed
fi
if [ ! -f gnomAD.del.bed ]
then
    zcat gnomad.v4.1.sv.sites.bed.gz | cut -f 1-5 | grep -w "DEL" | cut -f 1-4 | sort -k1,1V -k2,2n | uniq > gnomAD.del.bed
fi
if [ ! -f ont.vcf.gz ]
then
    bcftools view -O b -o tmp.bcf ${BASEDIR}/../release/svim-asm-hg38/svim.asm.hg38.bcf chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22
    bcftools index tmp.bcf
    bcftools annotate -c INFO -a ${BASEDIR}/../release/svim-asm-hg38/svim.asm.hg38.noGt.SVAN_1.3.bcf tmp.bcf | bgzip > ont.vcf.gz
    tabix ont.vcf.gz
    rm tmp.bcf tmp.bcf.csi
fi

WIN=100
echo -e "svclass\tsvtype\tdataset\ttotal\tshared\tunique\tnoveltyrate" > stats.tsv

## Insertions
for TAG in all
do
    echo ${TAG}
    echo "Insertions" `wc -l gnomad.ins.bed`
    bcftools view -i '(STRLEN(ALT)>=(STRLEN(REF)+50))' --min-ac 1 ont.vcf.gz | bcftools query -f "%CHROM\t%POS\t%END\t%ID\t%INFO/ITYPE_N\t%INFO/FAM_N\n" - | awk '{print $1"\t"($2-'${WIN}')"\t"($3+'${WIN}')"\t"$4"\t"$5"\t"$6;}' | sort -k1,1V -k2,2n | uniq > ont.${TAG}.bed
    COUNTONT=`cat ont.${TAG}.bed | cut -f 1-3 | sort | uniq | wc -l`
    SHAREDGAD=`bedtools intersect -a ont.${TAG}.bed -b <(zcat gnomad.v4.1.sv.sites.bed.gz  | cut -f 1-5 | grep -w "INS" | cut -f 1-4) -wao | awk '$7!="."' | cut -f 1-3 | sort | uniq | wc -l`
    UNIQGAD=`bedtools intersect -a ont.${TAG}.bed -b <(zcat gnomad.v4.1.sv.sites.bed.gz  | cut -f 1-5 | grep -w "INS" | cut -f 1-4) -wao | awk '$7=="."' | cut -f 1-3 | sort | uniq | wc -l`
    NRGAD=`echo "${UNIQGAD} / ${COUNTONT}" | bc -l`
    echo -e "${TAG}\tins\tgnomAD-SV\t${COUNTONT}\t${SHAREDGAD}\t${UNIQGAD}\t${NRGAD}" >> stats.tsv
    rm ont.${TAG}.bed
done

## MEIs
for TAG in Alu L1 SVA
do
    GTAG=`echo ${TAG} | sed 's/Alu/INS:ME:ALU/' | sed 's/L1/INS:ME:LINE1/' | sed 's/SVA/INS:ME:SVA/'`
    echo ${TAG} ${GTAG}
    echo "Insertions" ${TAG} `zcat gnomad.v4.1.sv.sites.bed.gz  | cut -f 1-5 | grep -w "${GTAG}" | cut -f 1-4 | sort | uniq | wc -l`
    bcftools view -i '(STRLEN(ALT)>=(STRLEN(REF)+50))' --min-ac 1 ont.vcf.gz | bcftools query -f "%CHROM\t%POS\t%END\t%ID\t%INFO/ITYPE_N\t%INFO/FAM_N\n" - | awk '{print $1"\t"($2-'${WIN}')"\t"($3+'${WIN}')"\t"$4"\t"$5"\t"$6;}' | egrep "partnered|solo" | grep -w "${TAG}" | sort -k1,1V -k2,2n | uniq > ont.${TAG}.bed
    COUNTONT=`cat ont.${TAG}.bed | cut -f 1-3 | sort | uniq | wc -l`
    SHAREDGAD=`bedtools intersect -a ont.${TAG}.bed -b <(zcat gnomad.v4.1.sv.sites.bed.gz  | cut -f 1-5 | grep -w "${GTAG}" | cut -f 1-4) -wao | awk '$7!="."' | cut -f 1-3 | sort | uniq | wc -l`
    UNIQGAD=`bedtools intersect -a ont.${TAG}.bed -b <(zcat gnomad.v4.1.sv.sites.bed.gz  | cut -f 1-5 | grep -w "${GTAG}" | cut -f 1-4) -wao | awk '$7=="."' | cut -f 1-3 | sort | uniq | wc -l`
    NRGAD=`echo "${UNIQGAD} / ${COUNTONT}" | bc -l`
    echo -e "${TAG}\tins\tgnomAD-SV\t${COUNTONT}\t${SHAREDGAD}\t${UNIQGAD}\t${NRGAD}" >> stats.tsv
    rm ont.${TAG}.bed
done

## Deletions
for TAG in all large
do
    echo "Deletions" `wc -l gnomAD.del.bed`
    bcftools view -i '(STRLEN(REF)>=(STRLEN(ALT)+50))' --min-ac 1 ont.vcf.gz | bcftools query -f "%CHROM\t%POS\t%END\t%ID\t%INFO/DTYPE_N\t%INFO/FAM_N\n" - | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6;}' | sort -k1,1V -k2,2n | uniq > ont.${TAG}.bed
    if [ ${TAG} == "large" ]
    then
	cat ont.${TAG}.bed | awk '($3-$2>250)' > ont.${TAG}.bed.tmp && mv ont.${TAG}.bed.tmp ont.${TAG}.bed
    fi
    COUNTONT=`cat ont.${TAG}.bed | cut -f 1-3 | sort | uniq | wc -l`
    SHAREDGAD=`python3 intersect.py -a ont.${TAG}.bed -b gnomAD.del.bed | cut -f 1-3 | sort | uniq | wc -l`
    UNIQGAD=`python3 intersect.py -v -a ont.${TAG}.bed -b gnomAD.del.bed | cut -f 1-3 | sort | uniq | wc -l`
    NRGAD=`echo "${UNIQGAD} / ${COUNTONT}" | bc -l`
    echo -e "${TAG}\tdel\tgnomAD-SV\t${COUNTONT}\t${SHAREDGAD}\t${UNIQGAD}\t${NRGAD}" >> stats.tsv
    rm ont.${TAG}.bed
done


## Left/Right extension
rm -f ext.tsv
for EXT in 1 10 50 100 200 500
do
    echo ${EXT}
    bcftools view -i '(STRLEN(ALT)>=(STRLEN(REF)+50))' --min-ac 1 ont.vcf.gz | bcftools query -f "%CHROM\t%POS\t%END\t%ID\t%INFO/ITYPE_N\t%INFO/FAM_N\n" - | awk '{print $1"\t"($2-'${EXT}')"\t"($3+'${EXT}')"\t"$4;}' | sort -k1,1V -k2,2n | uniq > ont.${EXT}.bed
    cat gnomad.ins.bed | awk '{print $1"\t"($2-'${EXT}')"\t"($3+'${EXT}')"\t"$4;}' > gnomad.${EXT}.bed
    COUNTONT=`cat ont.${EXT}.bed | cut -f 1-3 | sort | uniq | wc -l`
    COUNTAD=`cat gnomad.${EXT}.bed | cut -f 1-3 | sort | uniq | wc -l`
    SHAREDONT=`bedtools intersect -a ont.${EXT}.bed -b gnomad.${EXT}.bed -wao | awk '$5!="."' | cut -f 1-3 | sort | uniq | wc -l`
    UNIQONT=`expr ${COUNTONT} - ${SHAREDONT}`
    SHAREDAD=`bedtools intersect -a gnomad.${EXT}.bed -b ont.${EXT}.bed -wao | awk '$5!="."' | cut -f 1-3 | sort | uniq | wc -l`
    UNIQAD=`expr ${COUNTAD} - ${SHAREDAD}`
    NOVELONT=`echo "${UNIQONT} / ${COUNTONT}" | bc -l`
    echo -e "${EXT}\t1kG_ONT\t${COUNTONT}\t${SHAREDONT}\t${UNIQONT}\t${NOVELONT}" >> ext.tsv
    echo -e "${EXT}\tgnomAD\t${COUNTAD}\t${SHAREDAD}\t${UNIQAD}\tNA" >> ext.tsv
    rm ont.${EXT}.bed gnomad.${EXT}.bed
done


#Rscript novelty.R
