import sys
import argparse
from collections import defaultdict
import gzip

def parse_line(line, gt_list, pos, add_qual=None, pass_only=False):
	fields = line.split()
	var_id = '-'.join([fields[0], fields[1], fields[3], fields[4]])
	assert len(fields[4].split(',')) == 1
	assert var_id != '.'
	format_field = fields[8].split(':')
	genotype_field = fields[9].split(':')
	assert 'GT' in format_field
	index_gt = format_field.index('GT')
	gt = genotype_field[index_gt].replace('|', '/').split('/')
	assert gt[0] in ['.', '1', '0']
	gt_qual = 0
	qual_thresh = 0
	if add_qual != None:
		gt_qual = int(genotype_field[format_field.index('GQ')]) if genotype_field[format_field.index('GQ')] != '.' else 0
		qual_thresh = int(add_qual)
	if not '.' in gt:
		if len(gt) < 2:
			assert 'X' in fields[0] or 'Y' in fields[0]
			gt.append(gt[0])
		gt_list[var_id][pos] = int(gt[0]) + int(gt[1])
		gt_list[var_id][pos+3] = True
		if gt_qual < qual_thresh:
			gt_list[var_id][pos] = 3
		if pass_only and (fields[6] not in ['.', 'PASS']):
			gt_list[var_id][pos] = 0
			gt_list[var_id][pos+3] = False
	else:
		gt_list[var_id][pos] = 3
		gt_list[var_id][pos+3] = True


parser = argparse.ArgumentParser(prog='genotyping-evaluation.py', description=__doc__)
parser.add_argument('baseline', metavar='BASELINE', help='ground truth.')
parser.add_argument('callset', metavar='CALLSET', help='callset variants')
parser.add_argument('ids', metavar='IDS', help='genotyped ids')
parser.add_argument('--qual', metavar='QUALITY', default=None, help='quality threshold')
parser.add_argument('--only_pass', action='store_true', default=False, help="Only consider PASS or '.' calls.")
args = parser.parse_args()

total_baseline = 0
total_callset = 0
total_ids = 0
only_pass = args.only_pass
quality = int(args.qual) if args.qual else None
if quality == 0:
	quality = None

id_to_genotypes = defaultdict(lambda: [3, 0, None, False, False, False])
for id in open(args.ids, 'r'):
	id_to_genotypes[id.strip()] = [0, 0, None, False, False, True]
	total_ids += 1

for line in gzip.open(args.baseline, 'rt'):
	if line.startswith('#'):
		continue
	parse_line(line, id_to_genotypes, 0)
	total_baseline += 1

for line in gzip.open(args.callset, 'rt'):
	if line.startswith('#'):
		continue
	parse_line(line, id_to_genotypes, 1, add_qual=quality, pass_only=only_pass)
	total_callset += 1

sys.stderr.write('Read ' + str(total_baseline) + ' baseline variants.\n')
sys.stderr.write('Read ' + str(total_callset) + ' callset variants.\n') 
sys.stderr.write('Read ' + str(total_ids) + ' ids.\n')

intersections = {}
for a in [False, True]:
	for b in [False, True]:
		for c in [False, True]:
			intersections[(a,b,c)] = 0
print(intersections)

check_base = 0
check_call = 0
check_ids = 0

for var,stats in id_to_genotypes.items():
	intersections[(stats[3], stats[4], stats[5])] += 1

print('in_baseline\tin_callset\tin_ids')
for i,j in intersections.items():
	print('\t'.join([str(i[0]), str(i[1]), str(i[2]), str(j)]))
	if i[0]:
		check_base += j
	if i[1]:
		check_call += j
	if i[2]:
		check_ids += j

assert(check_base == total_baseline)
assert(check_call <= total_callset)
assert(check_ids == total_ids)

# construct matrix
matrix = [[0,0,0,0] for i in range(4)]

for var_id,genotypes in id_to_genotypes.items():
	matrix[genotypes[0]][genotypes[1]] += 1

print('\t'.join(['\t', '0', '1', '2', '3']))
for i,line in enumerate(matrix):
	print('\t'.join([str(i)] + [str(l) for l in line]))

all = sum(matrix[0]) + sum(matrix[1]) + sum(matrix[2])
all_present = sum(matrix[1]) + sum(matrix[2])
untyped = matrix[0][3] + matrix[1][3] + matrix[2][3]

correct_absent = matrix[0][0]
correct_het = matrix[1][1]
correct_hom = matrix[2][2]
total_absent = sum(matrix[0][0:3])
total_het = sum(matrix[1][0:3])
total_hom = sum(matrix[2][0:3])

conc_absent = correct_absent / float(total_absent) if total_absent > 0 else 0.0
conc_het = correct_het / float(total_het) if total_het > 0 else 0.0
conc_hom = correct_hom / float(total_hom) if total_hom > 0 else 0.0
balanced_concordance = (conc_absent + conc_het + conc_hom) / 3.0
overall_concordance = (correct_absent + correct_het + correct_hom) / (total_absent + total_het + total_hom) if (total_absent + total_het + total_hom) > 0 else 0.0
typed_all = 1.0 - (untyped / float(all)) if float(all) > 0 else 0.0

## vcfeval like stats for comparison
true_positives = correct_hom + correct_het
false_positives = matrix[0][1] + matrix[0][2] + matrix[1][2] + matrix[2][1] + matrix[3][1] + matrix[3][2]
false_negatives = matrix[1][0] + matrix[1][2] + matrix[1][3] + matrix[2][0] + matrix[2][1] + matrix[2][3]
precision = true_positives / float(true_positives + false_positives) if (true_positives + false_positives) > 0 else 0.0
recall = true_positives / float(true_positives + false_negatives) if (true_positives + false_negatives) > 0 else 0.0


unmatched = intersections[(False, True, False)]
print(unmatched, sum(matrix[3]))
assert sum(matrix[3]) >= unmatched


first_column = ['weighted_concordance', 'genotype_concordance_all', 'genotype_concordance_absent', 'genotype_concordance_het', 'genotype_concordance_hom', 'typed', 'total_baseline_variants', 'total_unmatched', 'precision', 'recall']
second_column = [str(balanced_concordance), str(overall_concordance), str(conc_absent), str(conc_het), str(conc_hom), str(typed_all), str(all_present), str(unmatched), str(precision), str(recall)]

for first,second in zip(first_column,second_column):
	print('\t'.join([first,second]))


