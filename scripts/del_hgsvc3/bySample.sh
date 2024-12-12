#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../mamba/bin:${PATH}

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
if [ ! -f trf.svs ]
then
    echo "TRF"
    zcat ont.vcf.gz | grep -v "^#" | awk '{print ">"$3"\n"$4;}' > ont.fasta
    ./trf409.linux64 ont.fasta 2 7 7 80 10 50 2000 -m -l 30 -h -ngs > ont.trf
    rm ont.fasta
    ./faCount ont.fasta.*.mask > ont.trf.faCount
    rm ont.fasta.*.mask ont.trf
    # Tandem repeat masked SVs (at least 50%)
    cat ont.trf.faCount | tail -n +2 | cut -f 1,2,7 | awk '$3/$2>0.5 {print $1;}' | sort | uniq > trf.svs
    rm ont.trf.faCount
fi
echo "Deletion type"
bcftools query -i '(STRLEN(REF)>=(STRLEN(ALT)+50))' -f "%INFO/DTYPE_N\n" ${BASEDIR}/../release/svan-annotation/final-vcf.unphased.SVAN_1.3.vcf.gz | sort | uniq -c

echo -e "size\tsample\tdataset\tall\tshared\tunique\tFDR\tsensitivity\tTRfraction" > stats.tsv
for SAMPLE in `cat samples.tsv`
do
    for TAG in small large
    do
	bcftools view -i '(STRLEN(REF)>=(STRLEN(ALT)+50))' -s ${SAMPLE} --min-ac 1 ont.vcf.gz | bcftools query -f "%CHROM\t%POS\t%END\t%ID\t%INFO/DTYPE_N\t%INFO/FAM_N\n" - | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6;}' | sort -k1,1V -k2,2n | uniq > ont.${SAMPLE}.${TAG}.bed
	if [ ${TAG} == "small" ]
	then
	    cat ont.${SAMPLE}.${TAG}.bed | awk '$3-$2<250' > ont.${SAMPLE}.${TAG}.bed.tmp && mv ont.${SAMPLE}.${TAG}.bed.tmp ont.${SAMPLE}.${TAG}.bed
	else
	    cat ont.${SAMPLE}.${TAG}.bed | awk '($3-$2>=250)' > ont.${SAMPLE}.${TAG}.bed.tmp && mv ont.${SAMPLE}.${TAG}.bed.tmp ont.${SAMPLE}.${TAG}.bed
	fi
	zcat variants_T2T-CHM13_sv_insdel_HGSVC2024v1.0.tsv.gz | grep -v "^chrX" | grep -v "^chrY" | cut -f 1-6,9,16,17 | grep "${SAMPLE}" | awk 'BEGIN {FS="\t"; IFS="\t";} (($5=="DEL") && ($6>=50)) {print $2"\t"$3"\t"$4"\t"$1;}' | sort -k1,1V -k2,2n | uniq > hgsvc.${SAMPLE}.${TAG}.bed
	COUNTONT=`cat ont.${SAMPLE}.${TAG}.bed | cut -f 1-3 | sort | uniq | wc -l`
	COUNTTRF=`cat ont.${SAMPLE}.${TAG}.bed | grep -w -Ff trf.svs | cut -f 1-3 | sort | uniq | wc -l`
	TRFS=`echo "${COUNTTRF} / ${COUNTONT}" | bc -l`
	COUNTHGSVC=`cat hgsvc.${SAMPLE}.${TAG}.bed | cut -f 1-3 | sort | uniq | wc -l`
	SHAREDONT=`python3 intersect.py -a ont.${SAMPLE}.${TAG}.bed -b hgsvc.${SAMPLE}.${TAG}.bed | cut -f 1-3 | sort | uniq | wc -l`
	UNIQONT=`expr ${COUNTONT} - ${SHAREDONT}`
	SHAREDHGSVC=`python3 intersect.py -a hgsvc.${SAMPLE}.${TAG}.bed -b ont.${SAMPLE}.${TAG}.bed | cut -f 1-3 | sort | uniq | wc -l`
	UNIQHGSVC=`expr ${COUNTHGSVC} - ${SHAREDHGSVC}`	
	echo ${SAMPLE} ${TAG}
	FDRONT=`echo "${UNIQONT} / ${COUNTONT}" | bc -l`
	SENSONT=`echo "${SHAREDHGSVC} / ${COUNTHGSVC}" | bc -l`
	FDRHGSVC=`echo "${UNIQHGSVC} / ${COUNTHGSVC}" | bc -l`
	SENSHGSVC=`echo "${SHAREDONT} / ${COUNTONT}" | bc -l`
	echo -e "${TAG}\t${SAMPLE}\t1kG_ONT\t${COUNTONT}\t${SHAREDONT}\t${UNIQONT}\t${FDRONT}\t${SENSONT}\t${TRFS}" >> stats.tsv
	echo -e "${TAG}\t${SAMPLE}\tHGSVC3\t${COUNTHGSVC}\t${SHAREDHGSVC}\t${UNIQHGSVC}\t${FDRHGSVC}\t${SENSHGSVC}\t${TRFS}" >> stats.tsv
	rm ont.${SAMPLE}.${TAG}.bed hgsvc.${SAMPLE}.${TAG}.bed 
    done
done

for TAG in small large
do
    # Avg. FDR
    echo "${TAG} average sites:" `cat stats.tsv  | grep "ONT" | grep "${TAG}" | awk '{SUM+=$4;} END {print SUM/NR;}'`
    echo "${TAG} average FDR:" `cat stats.tsv  | grep "ONT" | grep "${TAG}" | awk '{SUM+=$7;} END {print SUM/NR;}'`
    echo "${TAG} average TR fraction" `cat stats.tsv  | grep "ONT" | grep "${TAG}" | awk '{SUM+=$9;} END {print SUM/NR;}'`
done
