import argparse
import pandas
import sys
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
from scipy.stats import ttest_ind
from cyvcf2 import VCF
import copy

parser = argparse.ArgumentParser(prog='boxplot-lengths.py', description="Plot SV lengths from list of VCFs")
parser.add_argument('-meta', metavar='META', help='metadata')
parser.add_argument('-vcfs', metavar='VCFS', help='comma-separated list of vcfs')
parser.add_argument('-tables', metavar='TABLES', help='comma-separated tsv files generated using initial run with vcfs')
parser.add_argument('-names', metavar='NAMES', help='comma-separated list of names for the vcfs or tables')
parser.add_argument('-coverages', metavar='coverages', help='TSV file with coverages')
parser.add_argument('-output', metavar='OUTPUT', help='output directory')
args = parser.parse_args()

# reading metadata (sample-to-population map and color coding)
metadata = pandas.read_csv(args.meta, sep='\t', header=0)
metadata = metadata[["Sample name", "Population code", "Superpopulation code"]]

coverages = pandas.read_csv(args.coverages, sep='\t', header=0)
coverages = coverages[["SAMPLE", "T2T_median_cov"]]
sample2coverage = {}
for line in coverages.to_numpy():
    if line[1] > 15:
        sample2coverage[line[0]] = True
    else:
        sample2coverage[line[0]] = False

vcf_data = {}
names = args.names.split(',')

# reading vcfs and extracting genotypes
if args.vcfs:
    vcfs = args.vcfs.split(',')
    assert (len(names) == len(vcfs))
    for vcf, name in zip(vcfs, names):
        print("Reading", name)
        vcf = VCF(vcf)
        vcf_samples = vcf.samples
        data = {}
        for s in vcf_samples:
            sample_data=metadata[metadata["Sample name"] == s]
            data[s] = [[0,0,0], [0,0,0], sample_data["Population code"].values[0], sample_data["Superpopulation code"].values[0]]      # first element counts number of HET SVs [COMPLEX, DEL, INS]. second element counts HOM SVs [COMPLEX, DEL, INS].  
            
        sv_type_to_index = {'COMPLEX': 0, 'DEL': 1, 'INS': 2}
        # Writing data into file
        writer = open(args.output+'/boxplot-lengths-%s.tsv'%(name), 'w')
        print("Sample\tPopulation\tSuperpopulation\tHET_COMPLEX\tHET_DEL\tHET_INS\tHOM_COMPLEX\tHOM_DEL\tHOM_INS", file=writer)
        for variant in vcf:
            if variant.num_het+variant.num_hom_ref+variant.num_hom_alt == 0:
                continue
            _, _, sv_type, bub, length = variant.ID.split('-')
            length = int(length)
            if length < 50:
                continue
            for index, g in enumerate(variant.genotypes):
                sample = vcf_samples[index]
                if g[0] == -1:
                    continue
                #print(g[0], g[1])
                assert g[0] in [0, 1] and g[1] in [0, 1]
                if g[0] == 0 and g[1] == 0:
                    continue
                if g[0]+g[1] == 1:
                    # HET SV    
                    data[sample][0][sv_type_to_index[sv_type]] += length
                elif g[0]+g[1] == 2:
                    # HOM SV
                    data[sample][1][sv_type_to_index[sv_type]] += length
        vcf.close()
        for sample in vcf_samples:
            d = data[sample]
            print("%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d"%(sample, d[2], d[3], *d[0], *d[1]), file=writer)
        writer.close()
        vcf_data[name] = data
elif args.tables:
    tables = args.tables.split(',') 
    assert (len(names) == len(tables))
    for table, name in zip(tables, names):
        print("Reading", name)
        data = {}
        reader = open(table, 'r')
        for line in reader:
            if line[0] == "S":
                continue
            sample, pop, spop, het_com, het_del, het_ins, hom_com, hom_del, hom_ins = line.rstrip().split('\t')
            data[sample] = [[int(het_com), int(het_del), int(het_ins)], [int(hom_com), int(hom_del), int(hom_ins)], pop, spop]
        vcf_data[name] = data
else:
    sys.exit("No table or list of vcfs given")

# sort data by them superpopulation code
for name, data in vcf_data.items():
    vcf_data[name] = dict(sorted(data.items(), key=lambda x:x[1][3]))

color_by_pop = {'AFR': '#ffd845', 'AMR': '#710027', 'EAS': '#778500', 'EUR': '#018ead', 'SAS': '#c44cfd'}

# need to stratify data by population, coverage and variant type
het_data = {}
hom_data = {}
y = {'All': [], 'COMPLEX': [], 'DEL': [], 'INS': []}
pops = list(color_by_pop.keys())
covs = ['all', 'high', 'low']
for name in names:
    het_data[name] = {z: {x: copy.deepcopy(y) for x in pops} for z in covs}
    hom_data[name] = {z: {x: copy.deepcopy(y) for x in pops} for z in covs}

for name in names:
    het_data_by_name = het_data[name]
    hom_data_by_name = hom_data[name]
    data_by_name = vcf_data[name]
    for sample, value in data_by_name.items():
        cov = 'high' if sample2coverage[sample] else 'low'
        spop = value[3]
        for c in ['all', cov]:
            for data_dict, stored_data in zip([het_data_by_name, hom_data_by_name], [value[0], value[1]]):
                data_dict[c][spop]['All'].append(sum(stored_data)/1000000)
                data_dict[c][spop]['COMPLEX'].append(stored_data[0]/1000000)
                data_dict[c][spop]['DEL'].append(stored_data[1]/1000000)
                data_dict[c][spop]['INS'].append(stored_data[2]/1000000)

# create legend
handles = []
labels = []
n_names = len(names)
colors = sns.color_palette("tab10", n_names)

count = 0
for name in names:
    line = matplotlib.patches.Patch(color=colors[count], label=name)
    handles.append(line)
    labels.append(name)
    count += 1
count = 0

y_range = (0,7)
width = 0.9
tick_pos = [((n_names+1)*x)+(1+((n_names+1)/2)) for x in [0,1,2,3,4]]

print("Plotting")
for vtype in ['All', 'COMPLEX', 'DEL', 'INS']:
    for cov in covs:
        # HOM Counts for vtype and cov
        fig = plt.figure(figsize =(10, 10))
        bp = {}
        count = 0
        for name in names:
            data = [hom_data[name][cov][x][vtype] for x in pops]
            bp[name] = plt.boxplot(data, positions=[((n_names+1)*x)+(2+count) for x in [0,1,2,3,4]], patch_artist = True, notch=True, widths=width)
            for patch in bp[name]['boxes']:
                patch.set_facecolor(colors[count])
            for whisker in bp[name]['whiskers']:
                whisker.set(color ='black', linewidth = 1)
            for cap in bp[name]['caps']:
                cap.set(color ='black', linewidth = 1)
            for median in bp[name]['medians']:
                median.set(color ='black', linewidth = 1)
            for flier in bp[name]['fliers']:
                flier.set(marker ='o', color ='black', alpha = 1)
            count += 1
        plt.title("HOM LENGTH\ncoverage: %s   vtype: %s"%(cov, vtype))
        plt.xticks(tick_pos, pops, fontsize=15)
        plt.yticks(fontsize=15)
        plt.ylim(y_range)
        plt.xlim(1, ((n_names+1)*5)+1)
        plt.tight_layout()
        plt.legend(handles, labels)
        plt.savefig(args.output+'/boxplot-lengths-HOM-%s-%s.svg'%(cov, vtype))
        plt.close()

        # HET Counts for vtype and cov
        fig = plt.figure(figsize =(10, 10))
        bp = {}
        count = 0
        for name in names:
            data = [het_data[name][cov][x][vtype] for x in pops]
            bp[name] = plt.boxplot(data, positions=[((n_names+1)*x)+(2+count) for x in [0,1,2,3,4]], patch_artist = True, notch=True, widths=width)
            for patch in bp[name]['boxes']:
                patch.set_facecolor(colors[count])
            for whisker in bp[name]['whiskers']:
                whisker.set(color ='black', linewidth = 1)
            for cap in bp[name]['caps']:
                cap.set(color ='black', linewidth = 1)
            for median in bp[name]['medians']:
                median.set(color ='black', linewidth = 1)
            for flier in bp[name]['fliers']:
                flier.set(marker ='o', color ='black', alpha = 1)
            count += 1
        plt.title("HET LENGTH\ncoverage: %s   vtype: %s"%(cov, vtype))
        plt.xticks(tick_pos, pops, fontsize=15)
        plt.yticks(fontsize=15)
        plt.ylim(y_range)
        plt.xlim(1, ((n_names+1)*5)+1)
        plt.tight_layout()
        plt.legend(handles, labels)
        plt.savefig(args.output+'/boxplot-lengths-HET-%s-%s.svg'%(cov, vtype))
        plt.close()