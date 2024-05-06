import matplotlib.pyplot as plt
from cyvcf2 import VCF
import argparse
from collections import namedtuple, Counter
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

def process_SV(name):
    _, _, type, bub, len = name.split('-')
    return type, bub, len

def plot_rausch_curve(allele_counts, plot_title, output, fname, pp, subset):
    x = {}
    y = {}
    names = list(allele_counts.keys())
    n_names = len(names)
    
    for name in names:
        filtered_allele_counts = [a.ac for a in allele_counts[name][subset] if a.ac != 0]
        counter = Counter(filtered_allele_counts)
        sorted_counter = sorted(counter.items(), key=lambda pair: pair[0])
        x[name] = []
        y[name] = []
        for i,j in sorted_counter:
            x[name].append(i)
            y[name].append(j)
    
    fig, ax = plt.subplots(figsize = (10,10))
    count = 0
    colors = sns.color_palette("tab10", n_names)
    for name in names:
        plt.plot(x[name][:-1], y[name][:-1], color=colors[count], label=name)
        count += 1
    plt.xlabel("Variant Allele Count", fontsize=15)
    plt.ylabel("Number of Variant Sites", fontsize=15)
    plt.title(plot_title, fontsize=20)
    try:
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.xscale('log')
        plt.yscale('log')
        plt.tight_layout()
        plt.legend(fontsize=15)
        plt.savefig('%s-rausch-curve_%s.svg'%(output, fname))
        pp.savefig()
    except ValueError:
        pass
    plt.close()

SV = namedtuple('SV','name type bub len')
AC = namedtuple('AC','name type ac len')

parser = argparse.ArgumentParser(prog='rausch-curves.py', description="Plotting log-log rausch curves")
parser.add_argument('-vcfs', metavar='VCFS', help='comma-separated list of vcfs')
parser.add_argument('-names', metavar='NAMES', help='comma-separated list of names for vcfs')
parser.add_argument('-sample-sheet', metavar='sample_sheet', help='sample sheet containing list of unrelated samples')
parser.add_argument('-output', metavar='OUTPUT', help='output directory')
args = parser.parse_args()

# creating map of whether the sample is in the unrelated 908 list
unrelated_map = {}
sheet_reader = open(args.sample_sheet, 'r')
for line in sheet_reader:
    if line[0] == 's':
        continue
    sample, _, _, superpop, _, _, project = line.rstrip().split('\t')
    if project == '1kGP':
        unrelated_map[sample] = [True, superpop]
    else:
        unrelated_map[sample] = [False, superpop]

# initializing vcf
vcfs = args.vcfs.split(',')
names = args.names.split(',')
categories = ['all', 'unrelated', 'unrelated_AFR', 'unrelated_non-AFR', 'unrelated_AMR', 'unrelated_EAS', 'unrelated_EUR', 'unrelated_SAS']
data = {}
for vcf, name in zip(vcfs, names):
    vcf = VCF(vcf)
    print("Reading", name)
    samples = vcf.samples
    ac = {x: [] for x in categories}
    for variant in vcf:
        if variant.num_het+variant.num_hom_ref+variant.num_hom_alt == 0:
            continue
        id = variant.ID
        t, b, l = process_SV(id)
        if int(l) < 50:
            continue
        allele_count = {x: 0 for x in categories}
        for index, g in enumerate(variant.genotypes):
            if g[0] == -1:
                continue
            s = samples[index]
            allele_count['all'] += g[0] + g[1]
            if unrelated_map[s][0]:
                allele_count['unrelated'] += g[0] + g[1]
                allele_count['unrelated_%s'%(unrelated_map[s][1])] += g[0] + g[1]
                if unrelated_map[s][1] != 'AFR':
                    allele_count['unrelated_non-AFR'] += g[0] +g[1]
        for key in categories:
            ac[key].append(AC(id, t, allele_count[key], int(l)))
    vcf.close()
    data[name] = ac

print("Plotting")
pp = PdfPages('%s-rausch-curves.pdf'%(args.output))
# plotting the rausch curve for different variant sets
for subset in categories[1:2]:
    plot_rausch_curve(data, "%s Samples | All SVs with AC > 0"%(subset), args.output, '%s-SV-All'%(subset), pp, subset)
    '''
    for vtype in ['INS', 'DEL', 'COMPLEX']:
        new_data = {}
        for name in names:
            new_data[name] = {}
            new_data[name][subset] = [a for a in data[name][subset] if (a.type == vtype)]
        plot_rausch_curve(new_data, "%s Samples | All SV %s with AC > 0"%(subset, vtype), args.output, '%s-SV-%s'%(subset, vtype), pp, subset)
    '''
pp.close()