#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../mamba/bin:${PATH}

HG=${BASEDIR}/../genome/hg38.fa
export REF_CACHE=${BASEDIR}../genome/cache_hg38/%2s/%2s/%s

THREADS=8

for CRAM in ${BASEDIR}/../hg38/*.cram
do
    if [ -f ${CRAM} ]
    then
	ID=`echo ${CRAM} | sed 's/^.*\///' | sed 's/.hg38.cram$//'`
	echo ${ID}
	NanoPlot -t ${THREADS} --cram ${CRAM} -o ${ID}
    fi
done
