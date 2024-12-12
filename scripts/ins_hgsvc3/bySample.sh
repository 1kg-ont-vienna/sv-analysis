#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../mamba/bin:${PATH}

# Iterate Samples
WIN=100

if [ ! -f variants_T2T-CHM13_sv_insdel_HGSVC2024v1.0.tsv.gz ]
then
    wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/release/Variant_Calls/1.0/T2T-CHM13/annotation_table/variants_T2T-CHM13_sv_insdel_HGSVC2024v1.0.tsv.gz
fi

if [ ! -f ont.vcf.gz ]
then
    bcftools view -O b -o tmp.bcf -S samples.tsv ${BASEDIR}/../release/giggles-genotyping/final-vcf.unphased.vcf.gz chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22
    bcftools index tmp.bcf
    bcftools annotate -c INFO -a ${BASEDIR}/../release/svan-annotation/final-vcf.unphased.SVAN_1.3.vcf.gz tmp.bcf | bgzip > ont.vcf.gz
    tabix ont.vcf.gz
    rm tmp.bcf tmp.bcf.csi
fi
echo "All insertions >=50bp"
bcftools query -i '(STRLEN(ALT)>=(STRLEN(REF)+50))' -f "%INFO/ITYPE_N\t%INFO/FAM_N\n" ${BASEDIR}/../release/svan-annotation/final-vcf.unphased.SVAN_1.3.vcf.gz | wc -l
echo "Non tandem-repeat insertions"
bcftools query -i '(STRLEN(ALT)>=(STRLEN(REF)+50))' -f "%INFO/ITYPE_N\t%INFO/FAM_N\n" ${BASEDIR}/../release/svan-annotation/final-vcf.unphased.SVAN_1.3.vcf.gz | grep -v -P "^VNTR\t" | grep -v -P "^DUP\t" | wc -l

## ONT
echo -e "sample\tdataset\tall\tshared\tunique\tfdr\tsensitivity" > stats.tsv
for SAMPLE in `cat samples.tsv`
do
    # Insertions outside tandem repeats (covered by vamos)
    bcftools view -i '(STRLEN(ALT)>=(STRLEN(REF)+50)) && INFO/NOT_CANONICAL==0' -s ${SAMPLE} --min-ac 1 ont.vcf.gz | bcftools query -f "%CHROM\t%POS\t%ID\t%INFO/ITYPE_N\t%INFO/FAM_N\n" - | awk '{print $1"\t"$2"\t"($2+1)"\t"$3"\t"$4"\t"$5;}' | grep -v -P "\tVNTR\t" | grep -v -P "\tDUP\t" | sort -k1,1V -k2,2n | uniq > ont.${SAMPLE}.bed
    cat ont.${SAMPLE}.bed | cut -f 5 | sort | uniq -c

    # All insertions
    zcat variants_T2T-CHM13_sv_insdel_HGSVC2024v1.0.tsv.gz | grep -v "^chrX" | grep -v "^chrY" | cut -f 1-6,9,16,17 | grep "${SAMPLE}" | awk 'BEGIN {FS="\t"; IFS="\t";} (($5=="INS") && ($6>=50)) {print $2"\t"($3-'${WIN}')"\t"($4+'${WIN}')"\t"$1;}' | sort -k1,1V -k2,2n | uniq > hgsvc.${SAMPLE}.bed
    COUNTONT=`cat ont.${SAMPLE}.bed | cut -f 1-3 | sort | uniq | wc -l`
    COUNTHGSVC=`cat hgsvc.${SAMPLE}.bed | cut -f 1-3 | sort | uniq | wc -l`
    SHAREDONT=`bedtools intersect -a ont.${SAMPLE}.bed -b hgsvc.${SAMPLE}.bed -wao | awk '$7!="."' | cut -f 1-3 | sort | uniq | wc -l`
    UNIQONT=`bedtools intersect -a ont.${SAMPLE}.bed -b hgsvc.${SAMPLE}.bed -wao | awk '$7=="."' | cut -f 1-3 | sort | uniq | wc -l`
    SHAREDHGSVC=`bedtools intersect -a hgsvc.${SAMPLE}.bed -b ont.${SAMPLE}.bed -wao | awk '$5!="."' | cut -f 1-3 | sort | uniq | wc -l`
    UNIQHGSVC=`bedtools intersect -a hgsvc.${SAMPLE}.bed -b ont.${SAMPLE}.bed -wao | awk '$5=="."' | cut -f 1-3 | sort | uniq | wc -l`
    echo ${SAMPLE}
    FDRONT=`echo "${UNIQONT} / ${COUNTONT}" | bc -l`
    SENSONT=`echo "${SHAREDHGSVC} / ${COUNTHGSVC}" | bc -l`
    FDRHGSVC=`echo "${UNIQHGSVC} / ${COUNTHGSVC}" | bc -l`
    SENSHGSVC=`echo "${SHAREDONT} / ${COUNTONT}" | bc -l`
    echo -e "${SAMPLE}\t1kG_ONT\t${COUNTONT}\t${SHAREDONT}\t${UNIQONT}\t${FDRONT}\t${SENSONT}" >> stats.tsv
    echo -e "${SAMPLE}\tHGSVC3\t${COUNTHGSVC}\t${SHAREDHGSVC}\t${UNIQHGSVC}\t${FDRHGSVC}\t${SENSHGSVC}" >> stats.tsv
    rm ont.${SAMPLE}.bed hgsvc.${SAMPLE}.bed 
done

# Avg. FDR
echo "average FDR" `cat stats.tsv  | grep "ONT" | awk '{SUM+=$6;} END {print SUM/NR;}'`
