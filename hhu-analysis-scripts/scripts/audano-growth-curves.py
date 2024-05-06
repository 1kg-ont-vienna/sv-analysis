import random
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.legend_handler import HandlerTuple
import pandas
from cyvcf2 import VCF
import argparse
from collections import namedtuple, Counter

def process_SV(name):
    _, _, type, bub, len = name.split('-')
    return type, bub, len

def plot_audano_curves(svs_per_sample, map, metadata, output, prefix):
    # TODO: Create Legend. Figure out how to better represent the superpopulations when number of samples increases
    samples = list(svs_per_sample.keys())
    # create legend
    handles = []
    labels = []
    handles.append([mpatches.Patch(facecolor=c, label='Singletons (AC = 1)') for c in ['#ffcd33', '#bfbcb0']])
    labels.append('Singletons (AC = 1)')
    handles.append([mpatches.Patch(facecolor=c, label='Doubletons (AC = 2)') for c in ['#f5b900', '#96917d']])
    labels.append('Doubletons (AC = 2)')
    handles.append([mpatches.Patch(facecolor=c, label='Polymorphic SVs (AC > 2)') for c in ['#b88b00', '#96917d']])
    labels.append('Polymorphic SVs (AC > 2)')
    handles.append([mpatches.Patch(facecolor=c, label='Major SVs (AF >= 50%)') for c in ['#7a5c00', '#514e42']])
    labels.append('Major SVs (AF >= 50%)')
    handles.append(mpatches.Patch(facecolor='red', label='Shared SVs (AF = 100%)'))
    labels.append('Shared SVs (AF = 100%)')
    colors_1 = []
    colors_2 = []
    colors_3 = []
    colors_4 = []
    for s in samples:
        colors_1.append('#ffd845' if metadata[metadata['Sample name'] == s]['Superpopulation code'].values[0] == 'AFR' else '#bfbcb0')
        colors_2.append('#f5b900' if metadata[metadata['Sample name'] == s]['Superpopulation code'].values[0] == 'AFR' else '#9e9a87')
        colors_3.append('#b88b00' if metadata[metadata['Sample name'] == s]['Superpopulation code'].values[0] == 'AFR' else '#7a7563')
        colors_4.append('#7a5c00' if metadata[metadata['Sample name'] == s]['Superpopulation code'].values[0] == 'AFR' else '#514e42')
    count_shared_svs = []
    count_major_svs = []
    count_polymorphic_svs = []
    count_doubletons_svs = []
    count_singleton_svs = []
    sv_counter = Counter()
    
    count_shared_bubs = []
    count_major_bubs = []
    count_polymorphic_bubs = []
    count_doubletons_bubs = []
    count_singleton_bubs = []
    bub_counter = Counter()
    
    for n, sample in enumerate(samples):
        sample_svs = [sv.name for sv in svs_per_sample[sample] if int(sv.len) > 49]
        sample_bubbles = list(set([map[sv] for sv in sample_svs]))
        sv_counter.update(sample_svs)
        bub_counter.update(sample_bubbles)
        shared_counter = 0
        major_counter = 0
        polymorphic_counter = 0
        doubleton_counter = 0
        singleton_counter = 0
        for value in sv_counter.values():
            if value == n+1:
                singleton_counter += 1
                doubleton_counter += 1
                polymorphic_counter += 1
                major_counter += 1
                shared_counter += 1
            elif value >= (n+1)/2:
                singleton_counter += 1
                doubleton_counter += 1
                polymorphic_counter += 1
                major_counter += 1
            elif value > 2:
                singleton_counter += 1
                doubleton_counter += 1
                polymorphic_counter += 1
            elif value == 2:
                singleton_counter += 1
                doubleton_counter += 1
            elif value == 1:
                singleton_counter += 1
        count_shared_svs.append(shared_counter)
        count_major_svs.append(major_counter)
        count_polymorphic_svs.append(polymorphic_counter)
        count_doubletons_svs.append(doubleton_counter)
        count_singleton_svs.append(singleton_counter)
        
        shared_counter = 0
        major_counter = 0
        polymorphic_counter = 0
        doubleton_counter = 0
        singleton_counter = 0
        for value in bub_counter.values():
            if value == n+1:
                singleton_counter += 1
                doubleton_counter += 1
                polymorphic_counter += 1
                major_counter += 1
                shared_counter += 1
            elif value >= (n+1)/2:
                singleton_counter += 1
                doubleton_counter += 1
                polymorphic_counter += 1
                major_counter += 1
            elif value > 2:
                singleton_counter += 1
                doubleton_counter += 1
                polymorphic_counter += 1
            elif value == 2:
                singleton_counter += 1
                doubleton_counter += 1
            elif value == 1:
                singleton_counter += 1
        count_shared_bubs.append(shared_counter)
        count_major_bubs.append(major_counter)
        count_polymorphic_bubs.append(polymorphic_counter)
        count_doubletons_bubs.append(doubleton_counter)
        count_singleton_bubs.append(singleton_counter)
    
    fig, ax = plt.subplots(figsize = (10,5), dpi=200)
    ax.bar(samples, count_singleton_svs, color=colors_1)
    ax.bar(samples, count_doubletons_svs, color=colors_2)
    ax.bar(samples, count_polymorphic_svs, color=colors_3)
    ax.bar(samples, count_major_svs, color=colors_4)
    ax.bar(samples, count_shared_svs, color='red')
    ax.grid(axis='y', color='grey', linestyle='--', linewidth=0.5, alpha=0.5)
    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    plt.ylabel("Cumulative SVs", fontsize=15)
    plt.xlabel("Samples", fontsize=15)
    plt.yticks(fontsize=10)
    plt.legend(handles, labels, framealpha=1, frameon=True, loc='upper left', handler_map = {list: HandlerTuple(None)})
    plt.tight_layout()
    plt.savefig('%s/%ssv.svg'%(output, prefix))
    plt.savefig('%s/%ssv.pdf'%(output, prefix))
    plt.close()
    print('Stats for %s/%ssv.svg'%(output, prefix))
    print("Number of Shared SVs at the end of the curve: ", count_shared_svs[-1])
    print("Number of Major SVs at the end of the curve: ", count_major_svs[-1]-count_shared_svs[-1])
    print("Number of Polymorphic SVs at the end of the curve: ", count_polymorphic_svs[-1]-count_major_svs[-1])
    print("Number of Doubleton SVs at the end of the curve: ", count_doubletons_svs[-1]-count_polymorphic_svs[-1])
    print("Number of Singleton SVs at the end of the curve: ", count_singleton_svs[-1]-count_doubletons_svs[-1])

    fig, ax = plt.subplots(figsize = (10,5), dpi=200)
    ax.bar(samples, count_singleton_bubs, color=colors_1)
    ax.bar(samples, count_doubletons_bubs, color=colors_2)
    ax.bar(samples, count_polymorphic_bubs, color=colors_3)
    ax.bar(samples, count_major_bubs, color=colors_4)
    ax.bar(samples, count_shared_bubs, color='red')
    ax.grid(axis='y', color='grey', linestyle='--', linewidth=0.5, alpha=0.5)
    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    plt.ylabel("Cumulative Bubbles", fontsize=15)
    plt.xlabel("Samples", fontsize=15)
    plt.yticks(fontsize=10)
    plt.legend(handles, labels, framealpha=1, frameon=True, loc='upper left', handler_map = {list: HandlerTuple(None)})
    plt.tight_layout()
    plt.savefig('%s/%sbubbles.svg'%(output, prefix))
    plt.savefig('%s/%sbubbles.pdf'%(output, prefix))
    plt.close()
    print('Stats for %s/%sbubbles.svg'%(output, prefix))
    print("Number of Shared Bubs at the end of the curve: ", count_shared_bubs[-1])
    print("Number of Major Bubs at the end of the curve: ", count_major_bubs[-1]-count_shared_bubs[-1])
    print("Number of Polymorphic Bubs at the end of the curve: ", count_polymorphic_bubs[-1]-count_major_bubs[-1])
    print("Number of Doubleton Bubs at the end of the curve: ", count_doubletons_bubs[-1]-count_polymorphic_bubs[-1])
    print("Number of Singleton Bubs at the end of the curve: ", count_singleton_bubs[-1]-count_doubletons_bubs[-1])

SV = namedtuple('SV','name type bub len')
AC = namedtuple('AC','name type ac len')
SIZE_DIST = namedtuple('SIZED_DIST', 'type count len')

parser = argparse.ArgumentParser(prog='audano-growth-curves.py', description="Plotting growth curves")
parser.add_argument('-vcf', metavar='VCF', help='multisample biallelic vcf')
parser.add_argument('-meta', metavar='META', help='sample metadata')
parser.add_argument('-map', metavar='MAP', help='bubble ids to allele ids map')
parser.add_argument('-sample-sheet', metavar='sample_sheet', help='sample sheet containing list of unrelated samples')
parser.add_argument('-output', metavar='OUTPUT', help='output directory')
args = parser.parse_args()

# reading map information
bub_map = {}
map_reader = open(args.map, 'r')
for line in map_reader:
    if line[0] == "#":
        continue
    bub_id, allele_ids = line.rstrip().split('\t')
    for allele_id in allele_ids.split(','):
        bub_map[allele_id] = bub_id
map_reader.close()

# creating map of whether the sample is in the unrelated 908 list
unrelated_map = {}
sheet_reader = open(args.sample_sheet, 'r')
for line in sheet_reader:
    if line[0] == 's':
        continue
    sample, _, _, superpop, _, _, project = line.rstrip().split('\t')
    if project == '1kGP':
        unrelated_map[sample] = True
    else:
        unrelated_map[sample] = False

# initializing vcf
vcf = VCF(args.vcf)
vcf_samples = vcf.samples

# reading metadata (sample-to-population map and color coding)
metadata = pandas.read_csv(args.meta, sep='\t', header=0)
metadata = metadata[["Sample name", "Population code", "Superpopulation code"]]
metadata = metadata[metadata["Sample name"].isin(vcf_samples)].sort_values(by=["Superpopulation code"], ascending=False)

# sorted sample list required for Audano curve
non_AFR_samples = list(metadata[metadata["Superpopulation code"] != 'AFR']['Sample name'].values)
AFR_samples = list(metadata[metadata["Superpopulation code"] == 'AFR']['Sample name'].values)
random.shuffle(non_AFR_samples)
sorted_vcf_samples = non_AFR_samples + AFR_samples

# things to store
svs_per_sample = {}
svs_per_sample_unrelated = {}

for sample in sorted_vcf_samples:
    svs_per_sample[sample] = []
    if unrelated_map[sample]:
        svs_per_sample_unrelated[sample] = []

for variant in vcf:
    if variant.num_het+variant.num_hom_ref+variant.num_hom_alt == 0:
        continue
    name = variant.ID
    t, b, l = process_SV(name)
    sv = SV(name, t, b, int(l))
    allele_count = 0
    for index, g in enumerate(variant.genotypes):
        if g[0] == -1:
            continue    
        allele_count += g[0] + g[1]
        s = vcf_samples[index]
        if g[0]+g[1] >= 1:
            svs_per_sample[s].append(sv)
            if unrelated_map[s]:
                svs_per_sample_unrelated[s].append(sv)

# plotting the audano curve for the SVs
plot_audano_curves(svs_per_sample, bub_map, metadata, args.output, 'audano-all-')
plot_audano_curves(svs_per_sample_unrelated, bub_map, metadata, args.output, 'audano-unrelated-')
