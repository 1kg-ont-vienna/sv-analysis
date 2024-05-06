import sys
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
import seaborn as sns

outname = sys.argv[1]
threshold = int(sys.argv[2])

varlengths = []
total_alleles = []

for line in sys.stdin:
	if line.startswith('#'):
		continue
	fields = line.split()
	info_fields = { i.split('=')[0] : i.split('=')[1].strip() for i in fields[7].split(';') if '=' in i}
	alleles = [fields[3]] + [a for a in fields[4].split(',')]
	nr_alleles = len(alleles) - 1
	varlen = max([len(a) for a in alleles])
	varids = info_fields['ID'].split(',')
	assert len(varids) == (len(alleles) - 1)
	varlengths.append(varlen)
	total_alleles.append(nr_alleles)
	if (nr_alleles > threshold):
		start = int(fields[1])-1
		end = start + len(fields[3])
		print('\t'.join([fields[0], str(start), str(end)]))

ax = sns.jointplot(x=varlengths, y=total_alleles, kind='hex', bins='log', xscale = 'log')
ax.set_axis_labels("bubble length (= length of longest allele)", "number of alleles in bubble")
plt.colorbar()
plt.tight_layout()
plt.savefig(outname)
