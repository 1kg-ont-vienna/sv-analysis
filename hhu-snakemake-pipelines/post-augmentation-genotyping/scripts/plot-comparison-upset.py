import sys, argparse
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import upsetplot
import pandas as pd

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='plot_upset.py', description=__doc__)
    parser.add_argument('-t', '--table', required=True, help='Table to which annotations shall be added.')
    parser.add_argument('-o', '--output', required=True, help='output name prefix.')
    parser.add_argument('-n', '--names', nargs='+', default=[], required=True, help='Names of columns to include in plot.')
    parser.add_argument('-r', '--regions', nargs='+', default=[], help='For each region, plot variants in/outside that region.')
    args = parser.parse_args()
    
    df = pd.read_csv(args.table, sep='\t')
    df = df[df[args.names].any(axis=1)]
    

    numbers = defaultdict(lambda: defaultdict(lambda : 0))

    with PdfPages(args.output) as pdf:
        plt.figure()
        variants = df.groupby(by=args.names).size()
        fig = plt.figure()
        upsetplot.plot(variants, sort_by='cardinality', show_counts='%d')
        fig.tight_layout()
        pdf.savefig()
        plt.close()
        
        # one plot per variant type
        vartypes = ["INS", "DEL", "OTHER"]
        for vartype in vartypes:
            fig = plt.figure()
            df_sub = df[df["var_type"] == vartype]
            variants = df_sub.groupby(by=args.names).size()
            upsetplot.plot(variants, sort_by='cardinality', show_counts='%d')
            plt.suptitle(vartype)
            fig.tight_layout()
            pdf.savefig()
            plt.close()
    
            min_len = [0, 1000, 10000, 100000, 1000000]
            for i in range(len(min_len)-1):
                df_sub = df[(df["var_type"] == vartype) & (df["length"] >= min_len[i]) & (df["length"] < min_len[i+1])] 
                variants = df_sub.groupby(by=args.names).size()
                if df_sub.size == 0:
                    continue
                upsetplot.plot(variants, sort_by='cardinality', show_counts='%d')
                plt.suptitle(vartype + ' [' + str(min_len[i]) + ',' + str(min_len[i+1]) + ')')
                fig.tight_layout()
                pdf.savefig()
                plt.close()
    
        for index, row in df.iterrows():
            key = tuple([row[n] for n in args.names])
            numbers['all'][key] +=1

        for r in args.regions:
            plt.figure()
            df_sub = df[df[r]]
            variants = df_sub.groupby(by=args.names).size()
            upsetplot.plot(variants, sort_by='cardinality', show_counts='%d')
            plt.suptitle(r)
            pdf.savefig()
            plt.close()
            
            for index, row in df_sub.iterrows():
                key = tuple([row[n] for n in args.names])
                numbers[r][key] +=1
        
            plt.figure()
            df_sub = df[~df[r]]
            variants = df_sub.groupby(by=args.names).size()
            upsetplot.plot(variants, sort_by='cardinality', show_counts='%d')
            plt.suptitle('not ' + r)
            pdf.savefig()
            plt.close()

            for index, row in df_sub.iterrows():
                key = tuple([row[n] for n in args.names])
                numbers['not_' + r][key] +=1

    for r in args.regions:
        print('\t'.join(args.names + [r, 'not_' + r]))
        for key in numbers['all'].keys():
            set_elements = []
            for i,e in enumerate(key):
                if e:
                    set_elements.append('1')
                else:
                    set_elements.append('0')
            set_elements = ','.join(set_elements)
            percentage_in = numbers[r][key] / numbers['all'][key]
            percentage_out = numbers['not_' + r][key] / numbers['all'][key]
            print('\t'.join([set_elements] + [str(percentage_in), str(percentage_out)]))
