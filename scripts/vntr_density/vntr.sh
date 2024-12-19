#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../mamba/bin:${PATH}

if [ ! -f CHM13_notinalldifficultregions.bed.gz ]
then
    wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.5/CHM13@all/Union/CHM13_notinalldifficultregions.bed.gz
fi
if [ ! -f vamos.effMotifs-0.1.T2T-CHM13.tsv.gz ]
then
    wget -O vamos.effMotifs-0.1.T2T-CHM13.tsv.gz https://zenodo.org/records/13263615/files/vamos.effMotifs-0.1.T2T-CHM13.tsv.gz?download=1
fi

T2TSIZE=3093002071
echo -e "tandem_repeat_class\tgenomic_area\tarea_size\ttr_events\ttr_density"
for TR in VNTR STR
do
    BEDSIZE=`zcat CHM13_notinalldifficultregions.bed.gz | awk '{SUM+=($3-$2);} END {print SUM;}'`
    OV=`bedtools intersect -a <(zcat vamos.effMotifs-0.1.T2T-CHM13.tsv.gz  | awk '$6=="'${TR}'"' | cut -f 1-3) -b <(zcat CHM13_notinalldifficultregions.bed.gz) | sort | uniq | wc -l`
    DENS=`echo " ( ( ${OV} / ${BEDSIZE} ) * 1000000 ) " | bc -l`
    echo -e "${TR}\thigh-confident\t${BEDSIZE}\t${OV}\t${DENS}"
    BEDSIZE=`expr ${T2TSIZE} - ${BEDSIZE}`
    OV=`bedtools intersect -v -a <(zcat vamos.effMotifs-0.1.T2T-CHM13.tsv.gz  | awk '$6=="'${TR}'"' | cut -f 1-3) -b <(zcat CHM13_notinalldifficultregions.bed.gz) | sort | uniq | wc -l`
    DENS=`echo " ( ( ${OV} / ${BEDSIZE} ) * 1000000 ) " | bc -l`
    echo -e "${TR}\tdifficult\t${BEDSIZE}\t${OV}\t${DENS}"
done
