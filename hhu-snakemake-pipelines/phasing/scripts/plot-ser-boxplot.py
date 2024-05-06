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
parser.add_argument('-meta', metavar='META', help='sample metadata')
parser.add_argument('-output', metavar='OUTPUT', help='output directory')
args = parser.parse_args()

tsv_list = pandas.read_csv(args.tsvs, header=None).to_numpy()[:,0]
sample2pop, pop2color = read_metadata(args.meta)

sample_results = {}
for tsv in tsv_list:
    sample, data = read_tsv(tsv)
    in_trio = False if data[('trio', 'longread')][0] == 0 else True
    sample_results[sample] = [data, in_trio]

# these are the SER values for longread vs nygc phasing
non_trio_ser_values_by_pop = {'AFR': [], 'AMR': [], 'EAS': [], 'EUR': [], 'SAS': []}

for sample, data in sample_results.items():
    pop = sample2pop[sample]
    d = data[0][('longread', 'nygc')]
    ser = d[1]*100/d[0]
    non_trio_ser_values_by_pop[pop].append(ser)

## Some statistics
data = []
for pop in non_trio_ser_values_by_pop.keys():
    data += non_trio_ser_values_by_pop[pop]
data = numpy.array(data)
print("Mean SER percentage for all Longread vs NYGC comparison: ", data.mean())

## Plotting

# create legend
handles = []
labels = []
for pop,color in pop2color.items():
    line = matplotlib.patches.Patch(color=color, label=pop)
    handles.append(line)
    labels.append(pop)

fig = plt.figure(figsize =(10, 5))
ax = fig.add_subplot(111)
bp = ax.boxplot([non_trio_ser_values_by_pop[x] for x in non_trio_ser_values_by_pop.keys()], patch_artist = True)
colors = [v for _,v in pop2color.items()]
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)
for whisker in bp['whiskers']:
    whisker.set(color ='black', linewidth = 1)
for cap in bp['caps']:
    cap.set(color ='black', linewidth = 1)
for median in bp['medians']:
    median.set(color ='black', linewidth = 1)
for flier in bp['fliers']:
    flier.set(marker ='o', color ='black', alpha = 1)    
ax.set_xticklabels([x for x,_ in pop2color.items()])
plt.title("SER between NYGC Statistical Phasing and WhatsHap Longread Phasing", fontsize=15)
plt.ylabel("Switch Error Rate (%)", fontsize=15)
plt.yticks(fontsize=12)
plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
plt.tight_layout()
plt.legend(handles, labels, framealpha=1, frameon=True, fontsize=15)
plt.savefig(args.output+'/no-trio-ser.svg')
plt.savefig(args.output+'/no-trio-ser.pdf')