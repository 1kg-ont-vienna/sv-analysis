import sys
from variantclassifier import VariantType, determine_variant_from_line

import argparse

string_to_type = {
	"snp": [VariantType.snp],

	"small-insertion": [VariantType.small_insertion],
	"small-deletion": [VariantType.small_deletion],
	"small-complex": [VariantType.small_complex],

	"midsize-insertion": [VariantType.midsize_insertion],
	"midsize-deletion": [VariantType.midsize_deletion],
	"midsize-complex": [VariantType.midsize_complex],

	"large-insertion-1000": [VariantType.large_insertion_50_1000],
	"large-deletion-1000": [VariantType.large_deletion_50_1000],
	"large-complex-1000": [VariantType.large_complex_50_1000],

	"large-insertion-10000": [VariantType.large_insertion_1000_10000],
	"large-deletion-10000": [VariantType.large_deletion_1000_10000],
	"large-complex-10000": [VariantType.large_complex_1000_10000],

	"large-insertion-50000": [VariantType.large_insertion_10000_50000],
	"large-deletion-50000": [VariantType.large_deletion_10000_50000],
	"large-complex-50000": [VariantType.large_complex_10000_50000],

	"large-insertion-above50000": [VariantType.large_insertion_50000],
	"large-deletion-above50000": [VariantType.large_deletion_50000],
	"large-complex-above50000": [VariantType.large_complex_50000],

	"large-insertion": [VariantType.large_insertion_50_1000, VariantType.large_insertion_1000_10000, VariantType.large_insertion_10000_50000, VariantType.large_insertion_50000],
	"large-deletion": [VariantType.large_deletion_50_1000, VariantType.large_deletion_1000_10000, VariantType.large_deletion_10000_50000, VariantType.large_deletion_50000],
	"large-complex": [VariantType.large_complex_50_1000, VariantType.large_complex_1000_10000, VariantType.large_complex_10000_50000, VariantType.large_complex_50000],

	"indel": [VariantType.snp, VariantType.small_insertion, VariantType.small_deletion, VariantType.small_complex, VariantType.midsize_insertion, VariantType.midsize_deletion, VariantType.midsize_complex],
	"indels": [VariantType.small_insertion, VariantType.small_deletion, VariantType.small_complex, VariantType.midsize_insertion, VariantType.midsize_deletion, VariantType.midsize_complex],
	"sv": [VariantType.large_insertion_50_1000, VariantType.large_insertion_1000_10000, VariantType.large_insertion_10000_50000, VariantType.large_insertion_50000, VariantType.large_deletion_50_1000, VariantType.large_deletion_1000_10000, VariantType.large_deletion_10000_50000, VariantType.large_deletion_50000, VariantType.large_complex_50_1000, VariantType.large_complex_1000_10000, VariantType.large_complex_10000_50000, VariantType.large_complex_50000],
	
	"small": [VariantType.small_insertion, VariantType.small_deletion, VariantType.small_complex],
	"midsize" : [VariantType.midsize_insertion, VariantType.midsize_deletion, VariantType.midsize_complex],
	"all" : [VariantType.snp, VariantType.small_insertion, VariantType.small_deletion, VariantType.small_complex, VariantType.midsize_insertion, VariantType.midsize_deletion, VariantType.midsize_complex, VariantType.large_insertion_50_1000, VariantType.large_insertion_1000_10000, VariantType.large_insertion_10000_50000, VariantType.large_insertion_50000, VariantType.large_deletion_50_1000, VariantType.large_deletion_1000_10000, VariantType.large_deletion_10000_50000, VariantType.large_deletion_50000, VariantType.large_complex_50_1000, VariantType.large_complex_1000_10000, VariantType.large_complex_10000_50000, VariantType.large_complex_50000]
}

parser = argparse.ArgumentParser(prog='extract-varianttype.py', description=__doc__)
parser.add_argument('vartype', metavar='TYPE', help='variant type. Possible options: %s'%('|'.join(list(string_to_type.keys()))))
args = parser.parse_args()

assert args.vartype in string_to_type
type_to_extract = string_to_type[args.vartype]

total = 0
written = 0

for line in sys.stdin:
	if line.startswith('#'):
		print(line.strip())
		continue
	vartype = determine_variant_from_line(line)
	fields = line.strip().split()
	is_symbolic_vcf = fields[4].startswith('<')
	if vartype in type_to_extract or is_symbolic_vcf:
		print(line.strip())
		written += 1
	total += 1
sys.stderr.write('Extracted ' + str(written) + ' of ' + str(total) + ' variants.\n')
