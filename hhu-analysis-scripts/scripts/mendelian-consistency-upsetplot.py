import argparse
import upsetplot
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(prog='mendelian-consistency.py', description=__doc__)
parser.add_argument('-children', metavar='CHILDREN', required=True, help='Comma separated list of children.')
parser.add_argument('-out', metavar='OUT', required=True, help='output directory')

args = parser.parse_args()

children = args.children.split(',')

child_to_family_id = {'NA19828': '2418 (NA19828)',
                      'HG01258': 'CLM16 (HG01258)',
                      'HG00420': 'SH006 (HG00420)',
                      'NA19129': 'Y077 (NA19129)',
                      'NA12877': '1463 (NA12877)',
                      'NA12878': '1463 (NA12878)'}

inconsistent_variants = {}
for child in children:
    inconsistent_variants[child] = []
    bed_reader = open(args.out+'/'+child+'-all.bed', 'r')
    for record in bed_reader:
        _, _, _, id = record.rstrip().split('\t')
        inconsistent_variants[child].append(id)
    bed_reader.close()

upset_data = upsetplot.from_contents(inconsistent_variants)

fig = plt.figure(figsize=(20,10), dpi=200)
u = upsetplot.UpSet(upset_data, element_size=None)
u.plot(fig=fig)
plt.ylabel('# Variant Sites', fontsize=15)
plt.savefig(args.out+'/upsetplot-all-inconsistent_variant_sites.svg')
plt.savefig(args.out+'/upsetplot-all-inconsistent_variant_sites.pdf')
plt.close()

inconsistent_variants = {}
for child in children:
    inconsistent_variants[child] = []
    bed_reader = open(args.out+'/'+child+'-biallelic.bed', 'r')
    for record in bed_reader:
        _, _, _, id = record.rstrip().split('\t')
        inconsistent_variants[child].append(id)
    bed_reader.close()

upset_data = upsetplot.from_contents(inconsistent_variants)

fig = plt.figure(figsize=(20,10), dpi=200)
u = upsetplot.UpSet(upset_data, element_size=None)
u.plot(fig=fig)
plt.ylabel('# Variant Sites', fontsize=15)
plt.savefig(args.out+'/upsetplot-biallelic-inconsistent_variant_sites.svg')
plt.savefig(args.out+'/upsetplot-biallelic-inconsistent_variant_sites.pdf')
plt.close()

inconsistent_variants = {}
for child in children:
    inconsistent_variants[child] = []
    bed_reader = open(args.out+'/'+child+'-multiallelic.bed', 'r')
    for record in bed_reader:
        _, _, _, id = record.rstrip().split('\t')
        inconsistent_variants[child].append(id)
    bed_reader.close()

upset_data = upsetplot.from_contents(inconsistent_variants)

fig = plt.figure(figsize=(20,10), dpi=200)
u = upsetplot.UpSet(upset_data, element_size=None)
u.plot(fig=fig)
plt.ylabel('# Variant Sites', fontsize=15)
plt.savefig(args.out+'/upsetplot-multiallelic-inconsistent_variant_sites.svg')
plt.savefig(args.out+'/upsetplot-multiallelic-inconsistent_variant_sites.pdf')
plt.close()
