#!/bin/bash

cat summary.final-vcf.unphased_ALL.tsv | grep "FN$" | cut -f 1 | sed 's/-/\t/g'  | awk '{print $1"\t"$2"\t"($2+$4);}' | sed 's/^chr//' > in.table
Rscript fn.plot.R
rm in.table
