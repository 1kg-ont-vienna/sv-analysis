#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(prog='concordance-summary.py', description='summarizes concordance values into a tsv file.')
parser.add_argument('directory', metavar='DIRECTORY', help='directory to the concordance summaries.')
parser.add_argument('regions', metavar='REGIONS', help='comma-separated list of regions (can be multiallelic or biallelic).')
parser.add_argument('vartype', metavar='VARTYPE', help='comma-separated list of variant types.')
parser.add_argument('minaf', metavar='minaf', help='minimum allele frequency')
parser.add_argument('maxaf', metavar='maxaf', help='maximum allele frequency')

args = parser.parse_args()
regions = args.regions.split(',')
regions.sort()
vartype = args.vartype.split(',')
vartype.sort()
min_af = args.minaf
max_af = args.maxaf

title_row = ['variant_class','weighted_concordance', 'genotype_concordance_all', 'genotype_concordance_absent', 'genotype_concordance_het', 'genotype_concordance_hom', 'typed', 'total_baseline_variants', 'total_unmatched', 'precision', 'recall']
print('\t'.join(title_row))
for r in regions:
    for v in vartype:
        file = '%s/%s_%s_%s-%s/summary.txt'%(args.directory, r, v, min_af, max_af)
        result = ""
        result += "%s_%s"%(r,v)
        with open(file, 'r') as f:
            for n, line in enumerate(f):
                if n < 16:
                    continue
                val = float(line.rstrip().split('\t')[1])
                val = str(round(val, 4))
                result += "\t%s"%(val)
        print(result)
        
            

