# This uses the QC table created for the filtering step
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import sys
import pandas as pd
import argparse
from scipy.stats import pearsonr
from pylab import *
import seaborn as sns

def create_hwe_plot(df, outname, name):
    with PdfPages(outname + name + '.pdf') as pdf:
        for metric1, metric2 in zip(['PANEL_AF', 'CALLSET_UNRELATED_AF'], ['CALLSET_UNRELATED_AF', 'CALLSET_UNRELATED_HETEROZYGOSITY']):
            plt.figure()
            fig, ax = plt.subplots()
            x_values = []
            y_values = []
            for l,f in zip(df[metric1], df[metric2]):
                if not pd.isnull(f) and not pd.isnull(l):
                    x_values.append(l)
                    y_values.append(f)
            assert len(x_values) == len(y_values)
            if len(x_values) == 0:
                continue
            joint_kws=dict(gridsize=35, cmap="hot_r")
            ax = sns.jointplot(x=x_values, y=y_values, xlim=(-0.05, 1.05), ylim=(-0.05,1.05), bins='log', kind='hex', joint_kws=joint_kws, marginal_ticks=True, color="red")
            ax.set_axis_labels(metric1, metric2)
            if 'AF' in metric1 and 'AF' in metric2:
                pearson_corr, p_value = pearsonr(x_values, y_values)
                ax.fig.suptitle(name + " (r=" + str(pearson_corr) + ")")
            if  metric1 == 'CALLSET_UNRELATED_AF' and metric2 == 'CALLSET_UNRELATED_HETEROZYGOSITY':
                # plot theoretical line
                t = np.arange(0.0, 1.01, 0.01)
                s = [2*j*(1.0-j) for j in t]
                ax.fig.suptitle(name)
                ax.ax_joint.plot(t,s,'r-')
            plt.colorbar()
            plt.tight_layout()
            plt.savefig(outname + name + '-' + metric1 + '_' + metric2 + '.svg')
            pdf.savefig()
            plt.close()

parser = argparse.ArgumentParser(prog='hwe-af-simple.py', description="Making plots out of the QC table.")
parser.add_argument('-table', metavar='TABLE', help='statistics table')
parser.add_argument('-output', metavar='OUTPUT', help='output prefix')
args = parser.parse_args()

filename = args.table
outname = args.output

df = pd.read_csv(filename, sep='\t', header=0)

if 'CALLSET_UNRELATED_HETEROZYGOSITY' not in df.columns:
    df = df.assign(CALLSET_UNRELATED_HETEROZYGOSITY=lambda d: d['CALLSET_UNRELATED_NUM_HET_SAMPLES'] / d['CALLSET_UNRELATED_NUM_TOTAL_SAMPLES'])

svs = set(id for id in df.VARIANT_ID if (not 'SNV' in id) and (int(id.split('-')[-1])>=50))
large_insertions = set(id for id in svs if id.split('-')[2]=="INS")
large_deletions = set(id for id in svs if id.split('-')[2]=="DEL")
large_complex = set(id for id in svs if id.split('-')[2]=="COMPLEX")

for name, variants in zip(['large_deletions', 'large_insertions', 'large_complex', 'all'], [large_deletions, large_insertions, large_complex, svs]):
    df_sub=df[df.VARIANT_ID.isin(variants)].copy()
    df_sub.sort_values(by=['VARIANT_ID'], inplace=True)
    df_sub.set_index('VARIANT_ID')

    ids_autosomes = set([i for i in variants if not 'chrX' in i])

    df_autosomes = df_sub[df_sub.VARIANT_ID.isin(ids_autosomes)]
    
    create_hwe_plot(df_autosomes, outname=outname, name=name)