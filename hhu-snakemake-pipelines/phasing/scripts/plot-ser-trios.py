import argparse
import pandas
import numpy
import matplotlib.pyplot as plt
import matplotlib
from collections import namedtuple

def read_tsv(file):
    
    tsv = pandas.read_csv(file, sep='\t', header=0)
    sample = tsv['#sample'].to_numpy()[0]
    tsv = tsv[['chromosome', 'dataset_name0', 'dataset_name1', 'all_assessed_pairs', 'all_switches', 'blockwise_hamming']]
    data = {
        ('trio', 'longread'): numpy.array([0, 0, 0]),
        ('trio', 'trio-longread'): numpy.array([0, 0, 0]),
        ('trio', 'nygc'): numpy.array([0, 0, 0]),
        ('longread', 'trio-longread'): numpy.array([0, 0, 0]),
        ('longread', 'nygc'): numpy.array([0, 0, 0]),
        ('trio-longread', 'nygc'): numpy.array([0, 0, 0]),
        }
    
    for line in tsv.to_numpy():
        _, d1, d2, x, y, z = line
        data[(d1, d2)] += numpy.array([x, y, z])
    
    return sample, data

def read_metadata(file):
    metadata = pandas.read_csv(file, sep='\t', header=0)
    metadata = metadata[["Sample name", "Population code", "Superpopulation code"]]
    sample2pop = {}
    pop2color = {'AFR': '#FFCD33', 'AMR': '#FF3D3D', 'EAS': '#ADFF3D', 'EUR': '#64EBFF', 'SAS': '#FF30FF'}
    for line in metadata.to_numpy():
        sample, _, pop = line
        sample2pop[sample] = pop
    
    return sample2pop, pop2color

parser = argparse.ArgumentParser(prog='plot-ser-boxplot.py', description="Plotting the phasing comparison results")
parser.add_argument('-tsvs', metavar='TSVS', help='File containing list of Pairwise TSV output files from whatshap compare')
parser.add_argument('-output', metavar='OUTPUT', help='output directory')
args = parser.parse_args()

families = {
  '2418': 'NA19818_NA19819_NA19828'.split('_'),
  'CLM16': 'HG01256_HG01257_HG01258'.split('_'),
  'SH006': 'HG00418_HG00419_HG00420'.split('_'),
  'Y077': 'NA19128_NA19127_NA19129'.split('_'),
  '1463_Paternal': 'NA12889_NA12890_NA12877'.split('_'),
  '1463_Maternal': 'NA12891_NA12892_NA12878'.split('_')
}

sample2family = {}

for family, samples in families.items():
    for sample in samples:
        sample2family[sample] = family

tsv_list = pandas.read_csv(args.tsvs, header=None).to_numpy()[:,0]

sample_results = {}
for tsv in tsv_list:
    sample, data = read_tsv(tsv)
    sample_results[sample] = data

child_samples = ['NA19828', 'HG01258', 'HG00420', 'NA19129', 'NA12877', 'NA12878']

## Preparing data for plots

# switch error rate against the nygc phasing
longread_ser = {}
trio_ser = {}
triolongread_ser = {}

for sample in sample2family.keys():
    data = sample_results[sample]
    longread_ser[sample] = data[('longread', 'nygc')][1]*100/data[('longread', 'nygc')][0]
    trio_ser[sample] = data[('trio', 'nygc')][1]*100/data[('trio', 'nygc')][0]
    triolongread_ser[sample] = data[('trio-longread', 'nygc')][1]*100/data[('trio-longread', 'nygc')][0]

## Printing some stats

parent_data = []
child_data = []
for sample in longread_ser.keys():
        if sample in child_samples:
            child_data.append(longread_ser[sample])
        else:
            parent_data.append(longread_ser[sample])
print("Mean SER percentage for Longread phasing for parents: ", numpy.array(parent_data).mean())
print("Mean SER percentage for Longread phasing for child: ", numpy.array(child_data).mean())

parent_data = []
child_data = []
for sample in trio_ser.keys():
        if sample in child_samples:
            child_data.append(trio_ser[sample])
        else:
            parent_data.append(trio_ser[sample])
print("Mean SER percentage for Trio phasing for parents: ", numpy.array(parent_data).mean())
print("Mean SER percentage for Trio phasing for child: ", numpy.array(child_data).mean())

parent_data = []
child_data = []
for sample in triolongread_ser.keys():
        if sample in child_samples:
            child_data.append(triolongread_ser[sample])
        else:
            parent_data.append(triolongread_ser[sample])
print("Mean SER percentage for Trio+Longread phasing for parents: ", numpy.array(parent_data).mean())
print("Mean SER percentage for Trio+Longread phasing for child: ", numpy.array(child_data).mean())


## Plotting

fig = plt.figure(figsize =(10, 5))
x_labels = list(sample2family.keys())
x_pos = []
for i in range(len(families.keys())):
    x_pos.append(4*i+1)
    x_pos.append(4*i+2)
    x_pos.append(4*i+3)

plt.plot(x_pos, [y for y in longread_ser.values()], 'bo', label='Long-read phasing')
plt.plot(x_pos, [y for y in triolongread_ser.values()], 'go', label = 'Long-read + Trio phasing')
plt.plot(x_pos, [y for y in trio_ser.values()], 'ro', label = 'Trio phasing')
for x in [4*(i+1) for i in range(len(families.keys())-1)]:
    plt.axvline(x=x, color='black', ls='--', lw=1)
plt.xticks(x_pos, labels=x_labels, rotation=60, fontsize=12)
plt.yticks(fontsize=12)
plt.ylabel('Switch Error Rate (%)', fontsize=15)
plt.title('WhatsHap compare against NYGC Statistical Phasing', fontsize=15)
plt.legend(fontsize=15, framealpha=1)
plt.tight_layout()
plt.savefig(args.output+'/trio-ser.svg')
plt.savefig(args.output+'/trio-ser.pdf')