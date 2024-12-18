#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../mamba/bin:${PATH}
export SV_DIR=svtoolkit_1.04.1441/

rm -rf irsResults/
mkdir -p irsResults/

# Create IRS vcf
ID=DEL
BCF=${BASEDIR}/../release/shapeit5-phased-callset/final-vcf.phased.vcf.gz
if [ ! -f ${ID}.vcf ]
then
    bcftools view ${BCF} | head -n 500 | grep "^#" | sed 's/VCFv4.2/VCFv4.1/' > ${ID}.vcf
    bcftools view -i 'STRLEN(REF)>=(STRLEN(ALT)+50)' -m 2 -M 2 ${BCF} | grep -v "^#" | awk 'BEGIN {FS="\t"; OFS="\t";} {$6=100; $8="SVTYPE=DEL;END="($2+length($4)-length($5)); print $0;}' >> ${ID}.vcf
    echo ${ID} `cat ${ID}.vcf | grep -v "^#" | wc -l`
else
    # Run SAV
    java -Xmx8g -cp ${SV_DIR}/lib/SVToolkit.jar:${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar org.broadinstitute.sv.main.SVAnnotator -A IntensityRankSum -R /t2t.fa -vcf ${ID}.vcf -O irsResults/${ID}.affy.vcf -arrayIntensityFile ${BASEDIR}/../release/aCGH_1kGP_chm13/ALL.genome.Affy6_probe_intensity_matrix_2506samples.20120621.dat -sample ${BASEDIR}/../release/aCGH_1kGP_chm13/affy.sample.list -irsUseGenotypes true -writeReport True -reportFile irsResults/${ID}.affy.pval
fi
