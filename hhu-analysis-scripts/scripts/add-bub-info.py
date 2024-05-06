import argparse
import sys
import pandas
from cyvcf2 import VCF


parser = argparse.ArgumentParser(prog='add-bub-info.py', description="Adding bubble info for each biallelic variant record.")
parser.add_argument('-table', metavar='table', help='QC Table')
parser.add_argument('-panel', metavar='panel', help='Multiallelic panel VCF.')
parser.add_argument('-output', metavar='output', help='Output file.')
args = parser.parse_args()

table = pandas.read_csv(args.table, sep='\t', header=0)

panel = VCF(args.panel)

# mapping variant IDs to bubble IDs
var_id_to_bub = {}
bub_id_to_allele_len = {}
for variant in panel:
    bub_id = variant.ID
    var_ids = variant.INFO['ID']
    ref_length = len(variant.REF)
    alt_lengths = ','.join([str(len(x)) for x in variant.ALT])
    bub_id_to_allele_len[bub_id] = [ref_length, alt_lengths]
    for var_id in var_ids.split(','):
        for v in var_id.split(':'):
            if v in var_id_to_bub:
                assert bub_id == var_id_to_bub[v]
            var_id_to_bub[v] = bub_id

bub_id_list = []
ref_len_list = []
alt_len_list = []
for var_id in table['variant_id'].tolist():
    try:
        bub_id_list.append(var_id_to_bub[var_id])
        ref_len_list.append(bub_id_to_allele_len[var_id_to_bub[var_id]][0])
        alt_len_list.append(bub_id_to_allele_len[var_id_to_bub[var_id]][1])
    except KeyError:
        bub_id_list.append(None)
        ref_len_list.append(None)
        alt_len_list.append(None)
        print('No bubble id found for variant ', var_id, file=sys.stderr)

# renaming the columns which were incorrectly named in the previous script
table = table.rename(columns={"variant type": "variant_type", "variant length": "variant_length"})

table.insert(3, 'bub_id', bub_id_list)
table.insert(4, 'ref_allele_len', ref_len_list)
table.insert(5, 'alt_allele_len', alt_len_list)
table.to_csv(args.output, sep='\t', index=False)