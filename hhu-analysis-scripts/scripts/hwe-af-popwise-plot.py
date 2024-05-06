import seaborn
import matplotlib.pyplot as plt
import pandas
import argparse
import numpy as np
from scipy.stats import pearsonr

class Table:
  
    def __init__(self, stats) -> None:
        self._stats = stats
        pass
        
    def subset(self, var_type=None, is_SV=False, is_not_SV=False):
        assert not (is_SV == True and is_not_SV==True) 
        if var_type:
            assert var_type in ['SNV', 'INS', 'DEL', 'COMPLEX']
        if var_type:
            subtable = self._stats[self._stats['variant_type'] == var_type]
        else:
            subtable = self._stats
        if is_SV:
            subtable=subtable[subtable['variant_length'] >= 50]
        if is_not_SV:
            subtable=subtable[subtable['variant_length'] < 50]
        
        return subtable


def plot_all_hwe(data, prefix, out):
    '''
    Plots HWE curve for all the variant IDs for the entire cohort
    '''
    x_name=prefix+"_allele_freq"
    y_name=prefix+"_heterozygosity"
    title=out.split('/')[-1][0:-4]
    #g = seaborn.jointplot(data=data, x=x_name, y=y_name, kind='hex', height=10, ratio=5, ylim=(0, 1), cbar=True)
    
    g = seaborn.JointGrid(data=data, x=x_name, y=y_name, height=10, ratio=5, ylim=(0, 1))
    g.plot_marginals(seaborn.histplot, bins=50)
    g.plot_joint(plt.hexbin, bins='log', gridsize=50, cmap='Blues')
    g.set_axis_labels(x_name, y_name)
    x = np.linspace(0,1,10000)
    y = 2*x*(1-x)
    plt.plot(x,y, color='black', linewidth=2)
    plt.subplots_adjust(left=0, right=0.8, top=1, bottom=0)  # shrink fig so cbar is visible
    # make new ax object for the cbar
    plt.yticks(fontsize=15)
    plt.xticks(fontsize=15)
    cbar_ax = g.fig.add_axes([.85, .1, .05, .8])  # x, y, width, height
    plt.colorbar(cax=cbar_ax)
    plt.yticks(fontsize=15)
    g.fig.suptitle(title)
    g.savefig(out)
    plt.close()
    pass

def plot_pop_hwe(data, prefix, out):
    '''
    Plots HWE curve for all the variant IDs for the entire cohort
    '''
    x_name=prefix+"_allele_freq"
    y_name=prefix+"_heterozygosity"
    title=out.split('/')[-1][0:-4]
    g = seaborn.JointGrid(data=data, x=x_name, y=y_name, height=10, ratio=5, ylim=(0, 1))
    g.plot_marginals(seaborn.histplot, bins=50)
    g.plot_joint(plt.hexbin, bins='log', gridsize=50, cmap='Blues')
    g.set_axis_labels(x_name, y_name)
    x = np.linspace(0,1,10000)
    y = 2*x*(1-x)
    plt.plot(x,y, color='black', linewidth=2)
    plt.subplots_adjust(left=0, right=0.8, top=1, bottom=0)  # shrink fig so cbar is visible
    # make new ax object for the cbar
    plt.yticks(fontsize=15)
    plt.xticks(fontsize=15)
    cbar_ax = g.fig.add_axes([.85, .1, .05, .8])  # x, y, width, height
    plt.colorbar(cax=cbar_ax)
    plt.yticks(fontsize=15)
    g.fig.suptitle(title)
    g.savefig(out)
    plt.close()
    pass


def plot_all_panelAF_vs_callsetAF(data, prefix, out, writer):
    '''
    Plots panel allele freq vs callset allele freq for all the variant IDs
    '''
    x_name="panel_allele_freq"
    y_name=prefix+"_allele_freq"
    title=out.split('/')[-1][0:-4]
    print(title, file=writer)
    print(pearsonr(data[x_name], data[y_name]), file=writer)
    g = seaborn.JointGrid(data=data, x=x_name, y=y_name, height=10, ratio=5)
    g.plot_marginals(seaborn.histplot, bins=50)
    g.plot_joint(plt.hexbin, bins='log', gridsize=50, cmap='Blues')
    g.set_axis_labels(x_name, y_name)
    plt.subplots_adjust(left=0, right=0.8, top=1, bottom=0)  # shrink fig so cbar is visible
    # make new ax object for the cbar
    plt.yticks(fontsize=15)
    plt.xticks(fontsize=15)
    cbar_ax = g.fig.add_axes([.85, .1, .05, .8])  # x, y, width, height
    plt.colorbar(cax=cbar_ax)
    plt.yticks(fontsize=15)
    g.fig.suptitle(title)
    g.savefig(out)
    plt.close()
    pass

def plot_pop_panelAF_vs_callsetAF(data, prefix, out, writer):
    '''
    Plots panel allele freq vs callset allele freq for all the variant IDs
    '''
    x_name="panel_allele_freq"
    y_name=prefix+"_allele_freq"
    title=out.split('/')[-1][0:-4]
    print(title, file=writer)
    print(pearsonr(data[x_name], data[y_name]), file=writer)
    g = seaborn.JointGrid(data=data, x=x_name, y=y_name, height=10, ratio=5)
    g.plot_marginals(seaborn.histplot, bins=50)
    g.plot_joint(plt.hexbin, bins='log', gridsize=50, cmap='Blues')
    g.set_axis_labels(x_name, y_name)
    plt.subplots_adjust(left=0, right=0.8, top=1, bottom=0)  # shrink fig so cbar is visible
    # make new ax object for the cbar
    plt.yticks(fontsize=15)
    plt.xticks(fontsize=15)
    cbar_ax = g.fig.add_axes([.85, .1, .05, .8])  # x, y, width, height
    plt.colorbar(cax=cbar_ax)
    plt.yticks(fontsize=15)
    g.fig.suptitle(title)
    g.savefig(out)
    plt.close()
    pass

def print_basic_stats(table, writer):
    SVtable = table[table['variant_length'] >= 50]
    print("\tNumber of variants: ", table.shape[0], file=writer)
    print("\tNumber of SVs: ", table[table['variant_length'] >= 50].shape[0], file=writer)
    print("\tNumber of INS: ", table[table['variant_type'] == 'INS'].shape[0], file=writer)
    print("\tNumber of DEL: ", table[table['variant_type'] == 'DEL'].shape[0], file=writer)
    print("\tNumber of COMPLEX: ", table[table['variant_type'] == 'COMPLEX'].shape[0], file=writer)
    print("\tNumber of INS SV: ", SVtable[SVtable['variant_type'] == 'INS'].shape[0], file=writer)
    print("\tNumber of DEL SV: ", SVtable[SVtable['variant_type'] == 'DEL'].shape[0], file=writer)
    print("\tNumber of COMPLEX SV: ", SVtable[SVtable['variant_type'] == 'COMPLEX'].shape[0], file=writer)
    print("\tNumber of SNV: ", table[table['variant_type'] == 'SNV'].shape[0], file=writer)

def filtered_variant_stats(table: pandas.DataFrame, writer):
    bub_ids = set(table['bub_id'].tolist())
    print("\tNumber of variants filtered out: ", table.shape[0], file=writer)
    print("\tNumber of bubbles involved: ", len(bub_ids), file=writer)
    alt_lengths = table[['bub_id', 'alt_allele_len']]
    alt_lengths = alt_lengths.drop_duplicates()
    n_alleles = []
    for alt in alt_lengths['alt_allele_len'].tolist():
        n_alleles.append(len(alt.split(',')))
    
    print(pandas.Series(n_alleles).describe(), file=writer)

parser = argparse.ArgumentParser(prog='hwe-af-popwise-plots.py', description="Making plots out of the VCF callset stats population and variant type stratified.")
parser.add_argument('-table', metavar='TABLE', help='statistics table')
parser.add_argument('-output', metavar='OUTPUT', help='output prefix')
args = parser.parse_args()

writer = open(args.output+'.stdout', 'w')

table = pandas.read_csv(args.table, sep='\t', header=0)

table = table[table['allpop-all_allele_freq'] != 0]

# creating Table object for plotting
stats = Table(table)

# making a dataframe with only SVs
stats_SV = stats.subset(is_SV=True)
# plotting HWE for all samples for SVs only
plot_all_hwe(stats_SV, "allpop-all", args.output+'-all-hwe.SVonly.png')
# plotting AF plot for all samples for all variant types
plot_all_panelAF_vs_callsetAF(stats_SV, "allpop-unrelated", args.output+'-all-af.SVonly.png', writer)

# plotting HWE for unrelated samples for SVs only
plot_all_hwe(stats_SV, "allpop-all", args.output+'-unrelated-hwe.SVonly.png')
# plotting AF plot for unrelated samples for all variant types
plot_all_panelAF_vs_callsetAF(stats_SV, "allpop-unrelated", args.output+'-unrelated-af.SVonly.png', writer)


# plotting HWE and AF for pop-wise samples for SVs only
for pop in ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']:
    for u in ['all', 'unrelated']:
        plot_pop_hwe(stats_SV, pop+"-"+u, args.output+'-%s-hwe.%s.SVonly.png'%(u, pop))
        plot_pop_panelAF_vs_callsetAF(stats_SV, pop+"-"+u, args.output+'-%s-af.%s.SVonly.png'%(u, pop), writer)



# SV only
stats_SV_vtype = {}
stats_SV_vtype['INS'] = stats.subset(is_SV=True, var_type='INS')
stats_SV_vtype['DEL'] = stats.subset(is_SV=True, var_type='DEL')
stats_SV_vtype['COMPLEX'] = stats.subset(is_SV=True, var_type='COMPLEX')

# plotting HWE and AF
for vtype in ['INS', 'DEL', 'COMPLEX']:
    plot_all_hwe(stats_SV_vtype[vtype], "allpop-all", args.output+'-all-hwe.%s.SVonly.png'%(vtype))
    plot_all_panelAF_vs_callsetAF(stats_SV_vtype[vtype], "allpop-all", args.output+'-all-af.%s.SVonly.png'%(vtype), writer)
    plot_all_hwe(stats_SV_vtype[vtype], "allpop-unrelated", args.output+'-unrelated-hwe.%s.SVonly.png'%(vtype))
    plot_all_panelAF_vs_callsetAF(stats_SV_vtype[vtype], "allpop-unrelated", args.output+'-unrelated-af.%s.SVonly.png'%(vtype), writer)

writer.close()