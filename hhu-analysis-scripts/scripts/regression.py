import sys
import pandas as pd
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
import upsetplot
from scipy.stats import pearsonr
from pylab import *
from collections import defaultdict

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.ensemble import AdaBoostRegressor
from sklearn.svm import SVR
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer



if __name__ == "__main__":
    filename = sys.argv[1]
    outname = sys.argv[2]

    df = pd.read_csv(filename, sep='\t', header=0)
    df = df.assign(mendelian_consistency=lambda d: d['MENDEL_INCONSISTENT_TRIO'] / d['MENDEL_TOTAL_UNTRIVIAL_TRIO'])
    df = df.assign(CALLSET_HETEROZYGOSITY=lambda d: d['CALLSET_NUM_HET_SAMPLES'] / d['CALLSET_NUM_TOTAL_SAMPLES'])
    df = df.assign(CALLSET_UNRELATED_HETEROZYGOSITY=lambda d: d['CALLSET_UNRELATED_NUM_HET_SAMPLES'] / d['CALLSET_UNRELATED_NUM_TOTAL_SAMPLES'])

    df = df.assign(ac0_fail = lambda df: ~(df['CALLSET_AC'] > 0))
    df = df.assign(mendel_fail = lambda df: ((df['MENDEL_INCONSISTENT_TRIO'] >= 1)))
    df = df.assign(maxr2_fail = lambda df: df['PHASED_MAXR2'] <= 0.5)
    df = df.assign(hwe_fail = lambda df: (df['HWE_CHI_P_VALUE'] < 0.05))
    df = df.assign(carrier_gq_fail= lambda df: df.CALLSET_UNRELATED_CARRIER_AVG_GQ<50)
    df = df.assign(non_carrier_gq_fail= lambda df: df.CALLSET_UNRELATED_NON_CARRIER_AVG_GQ<50)
    
    df = df.assign(all_pass = lambda df: (df.PHASED_MAXR2 >= 0.8) & (df.MENDEL_UNTRIVIAL_CONSISTENT_TRIO >= 1) & (df.HWE_CHI_P_VALUE>0.05))
    df = df.assign(negative = lambda df: (((1*df.mendel_fail) + (1*df.maxr2_fail) + (1*df.hwe_fail) + (1*df.carrier_gq_fail) + (1*df.non_carrier_gq_fail)) >= 3))
 
    svs = set(id for id in df.VARIANT_ID if (not 'SNV' in id) and (int(id.split('-')[-1])>=50))

    large_insertions = set(id for id in svs if id.split('-')[2]=="INS")
    large_deletions = set(id for id in svs if id.split('-')[2]=="DEL")
    large_complex = set(id for id in svs if id.split('-')[2]=="COMPLEX")
    assert len(large_insertions) + len(large_deletions) + len(large_complex) == len(svs)

    # for pangenie, consider allele frequencies computed across unrelated samples (+ panel samples) only, and skip children (unless they are panel samples)
    metrics = ['PANEL_AF', 'CALLSET_UNRELATED_AF', 'CALLSET_UNRELATED_HETEROZYGOSITY']

    # features used for regression
    # not included anything from mendelian or self genotyping (but they are used for the strict set)
    features = [
        'MULTIALLELIC',
        'PANEL_AF',
        'PANEL_AN',
        'PANEL_AC_PSEUDO',
        'CALLSET_AF',
        'CALLSET_AN',
        'CALLSET_NUM_UNTYPED_GENOTYPES',
        'CALLSET_NUM_SAMPLES_GQ>=20',
        'CALLSET_NUM_SAMPLES_GQ>=50',
        'CALLSET_NUM_SAMPLES_GQ>=100',
        'CALLSET_NUM_SAMPLES_GQ>=250',
        'CALLSET_HETEROZYGOSITY',
        'CALLSET_CARRIER_AVG_GQ',
        'CALLSET_NON_CARRIER_AVG_GQ',
        'HWE_EXACTTEST_STATISTIC',
        'HWE_CHI_TEST_STAT',
        'HWE_CHI_P_VALUE',
        'PHASED_AC',
        'PHASED_HWE',
        'PHASED_EXCHET',
        'PHASED_IRS',
        'PHASED_MAXR2'
    ]

    variant_types = {}
    # only consider non-empty sets of variants for analysis
    for name, variants in zip(['large_deletions', 'large_insertions', 'large_complex'], [large_deletions, large_insertions, large_complex]):
        if variants:
            variant_types[name] = [variants, ['unfiltered', 'lenient_0.0', 'lenient_-0.5', 'lenient_0.5', 'strict'], ['all-regions']]
    print('Considering variant classes: ' + ','.join([v for v in variant_types.keys()]))	

    ### Writing notes
    note_writer = open(outname+'_notes.txt', 'w')
    print("Regression type: regressor = SVR(C=50)", file=note_writer)
    print("\nFeatures for Regression:", file=note_writer)
    for f in features:
        print("\t"+f, file=note_writer)
    print("\nMetrics:", file=note_writer)
    for m in metrics:
        print("\t"+m, file=note_writer)
    print("\nFilters:", file=note_writer)
    print('''df = df.assign(ac0_fail = lambda df: ~(df['CALLSET_AC'] > 0))
    df = df.assign(mendel_fail = lambda df: ((df['MENDEL_INCONSISTENT_TRIO'] >= 1)))
    df = df.assign(maxr2_fail = lambda df: df['PHASED_MAXR2'] <= 0.5)
    df = df.assign(hwe_fail = lambda df: (df['HWE_CHI_P_VALUE'] < 0.05))
    df = df.assign(carrier_gq_fail= lambda df: df.CALLSET_UNRELATED_CARRIER_AVG_GQ<50)
    df = df.assign(non_carrier_gq_fail= lambda df: df.CALLSET_UNRELATED_NON_CARRIER_AVG_GQ<50)
    
    df = df.assign(all_pass = lambda df: (df.PHASED_MAXR2 >= 0.8) & (df.MENDEL_UNTRIVIAL_CONSISTENT_TRIO >= 1) & (df.HWE_CHI_P_VALUE>0.05))
    df = df.assign(negative = lambda df: (((1*df.mendel_fail) + (1*df.maxr2_fail) + (1*df.hwe_fail) + (1*df.carrier_gq_fail) + (1*df.non_carrier_gq_fail)) >= 3))''', file=note_writer)
    print("\nNumber of variants in positive set: ", df[df['all_pass'] == True].shape[0], file=note_writer)
    print("Number of variants in negative set: ", df[df['negative'] == True].shape[0], file=note_writer)

    note_writer.close()
    ###

    #exit()
    
    df['score_SVR'] = np.nan
    for name, info in variant_types.items():
        if 'large' in name:
            df_sub = df[df.VARIANT_ID.isin(info[0]) & ~df.ac0_fail].copy()
            df_sub.sort_values(by=['VARIANT_ID'], inplace=True)
            df_sub.set_index('VARIANT_ID')

            # impute missing values
            imp = IterativeImputer(max_iter=10, verbose=0, random_state=0)
            imp.fit(df_sub[features])
            imputed_df = pd.DataFrame(imp.transform(df_sub[features]), columns=df_sub[features].columns, index=df_sub.index)
            df_sub.update(imputed_df, overwrite=True)

            # Fit transform using all data points (labeled + unlabeled)
            df_sub = df_sub.assign(target = lambda df: (-1*df.negative) + (1*df.all_pass))
            x = df_sub.loc[:,features].values
            scaler = StandardScaler()
            scaler.fit_transform(x)

            # Train model only on labeled points
            autosomal_ids = [i for i in info[0] if not 'chrX' in i]
            df_labeled = df_sub[df.VARIANT_ID.isin(autosomal_ids) & (df_sub.target != 0)]
            x = scaler.transform(df_labeled.loc[:,features].values)
            y = df_labeled.loc[:,['target']].values
            print('Training regression model')
            regressor = SVR(C=50)
            #regressor = RandomForestRegressor(n_estimators=10, random_state=0)
            #regressor = AdaBoostRegressor(n_estimators=1, loss='square', random_state=0, learning_rate=0.1)
            #regressor = GradientBoostingRegressor(n_estimators=50, random_state=0, validation_fraction=0.2, n_iter_no_change=5)
            regressor.fit(x,y.ravel())

            # Apply to unlabeled data
            y_pred = regressor.predict(scaler.transform(df_sub.loc[:,features].values))
            column_label = "score_" + name
            df_score = pd.DataFrame({"VARIANT_ID":df_sub.VARIANT_ID, column_label:y_pred})

            # Add column with variant specific scores to table
            df = df.merge(df_score, on='VARIANT_ID', how='left')

        for filter in info[1]:
            variant = info[0]
            for region in info[2]:
                # collect all variants in the region
                if region == "all-regions":
                    ids_region = variant
                
                # create plots
                with PdfPages(outname + '_' + name + '_' + filter + '_' + region + '.pdf') as pdf:
                    ids_autosomes = set([i for i in ids_region if not 'chrX' in i])
                
                    if filter == 'unfiltered':
                        df_sub = df[df.VARIANT_ID.isin(ids_region)]
                        df_sub_autosomes = df[df.VARIANT_ID.isin(ids_autosomes)]
                    elif filter.startswith('lenient'):
                        for column in ['score_'+name, 'VARIANT_LENGTH', 'PANEL_AF', 'PANEL_AC', 'PANEL_AC_PSEUDO', 'CALLSET_AC', 'CALLSET_AF', 'CALLSET_AN', 'CALLSET_NUM_UNTYPED_GENOTYPES', 'CALLSET_AVG_GQ', 'CALLSET_CARRIER_AVG_GQ', 'CALLSET_NON_CARRIER_AVG_GQ', 'HWE_EXACTTEST_STATISTIC', 'HWE_CHI_TEST_STAT', 'HWE_CHI_P_VALUE', 'PHASED_HWE', 'PHASED_EXCHET', 'MULTIALLELIC', 'PHASED_IRS', 'PHASED_MAXR2']:
                            plt.figure()
                            fig, ax = plt.subplots()
                            df[df.VARIANT_ID.isin(variant) & df.negative].hist(ax=ax, column=column, bins=64, bottom = 0.1, alpha=0.5, color='red')
                            df[df.VARIANT_ID.isin(variant) & df.all_pass].hist(ax=ax, column=column, bins=64, bottom = 0.1, alpha=0.5, color='blue')
                            df[df.VARIANT_ID.isin(variant) & ~df.negative & ~df.all_pass & ~df.ac0_fail].hist(ax=ax, column=column, bins=64, bottom = 0.1, alpha=0.5, color='grey')
                            ax.set_yscale('log')
                            pdf.savefig()
                            plt.close()
    
                        score_cutoff = float(filter.split('_')[1])
                        df_sub = df[df.VARIANT_ID.isin(ids_region) & ( df.all_pass | ((~df.ac0_fail) & (df['score_'+name] > score_cutoff)) ) ]
                        df_sub_autosomes = df[df.VARIANT_ID.isin(ids_autosomes) & ( df.all_pass | ((~df.ac0_fail) & (df['score_'+name] > score_cutoff)) ) ]
                    else:
                        assert(filter == 'strict')
                        df_sub = df[df.VARIANT_ID.isin(ids_region) & df.all_pass]
                        df_sub_autosomes = df[df.VARIANT_ID.isin(ids_autosomes) & df.all_pass]

                    print('  variant count ' + filter + ' ' + region + ' ' + name + ':', len(df_sub))

                    """ # create upset plot with filters
                    filter_counts = df_sub.groupby(by=['ac0_fail','mendel_fail','gq_fail', 'self_fail']).size()
                    plt.figure()
                    upsetplot.plot(filter_counts, sort_by='cardinality')
                    pdf.savefig()
                    plt.close() """

                    # heatmaps
                    for i in range(len(metrics)):
                        for j in range(i+1, len(metrics)):
                                plt.figure()
                                fig, ax = plt.subplots()
                                x_values = []
                                y_values = []
                                for l,f in zip(df_sub_autosomes[metrics[i]], df_sub_autosomes[metrics[j]]):
                                    if not pd.isnull(f) and not pd.isnull(l):
                                        x_values.append(l)
                                        y_values.append(f)
                                assert len(x_values) == len(y_values)
                                if len(x_values) == 0:
                                    continue
    #							cmap = cm.get_cmap('Greys', 6) 
    #							joint_kws=dict(gridsize=35, cmap=cmap)
                                joint_kws=dict(gridsize=35, cmap="hot_r")
                                ax = sns.jointplot(x=x_values, y=y_values, xlim=(-0.05, 1.05), ylim=(-0.05,1.05), bins='log', kind='hex', joint_kws=joint_kws, marginal_ticks=True, color="red")
                                ax.set_axis_labels(metrics[i], metrics[j])
                                if 'AF' in metrics[i] and 'AF' in metrics[j]:
                                    pearson_corr, p_value = pearsonr(x_values, y_values)
                                    ax.fig.suptitle(name + " (r=" + str(pearson_corr) + ")")
                                    print('  pearson correlation ' + name + ' ' + filter + ' ' + region + ':', metrics[i], metrics[j], pearson_corr)
                                if 'CALLSET_UNRELATED_AF' in [metrics[i], metrics[j]] and 'CALLSET_UNRELATED_HETEROZYGOSITY' in [metrics[i], metrics[j]]:
                                    # plot theoretical line
                                    t = np.arange(0.0, 1.01, 0.01)
                                    s = [2*j*(1.0-j) for j in t]
                                    ax.fig.suptitle(name)
                                    ax.ax_joint.plot(t,s,'r-')
                                plt.colorbar()
                                plt.tight_layout()
                                pdf.savefig()
                                plt.close()
            
                    # boxplots
                    pop_metrics = ['CALLSET_UNRELATED_AF', 'PANEL_AF']
                    length_intervals = []
                    n_intervals = 10
                    previous = -0.1
                    length_intervals = []
                    for i in range(1, n_intervals + 1):
                        length_intervals.append([previous, i/n_intervals + 0.0001])
                        previous = i/n_intervals
                    print(length_intervals)
                    length_afs = [[] for i in range(n_intervals)]
                    for l,f in zip(df_sub_autosomes[pop_metrics[0]], df_sub_autosomes[pop_metrics[1]]):
                        if not pd.isnull(l):
                            for i,interval in enumerate(length_intervals):
                                if interval[0] < f <= interval[1]:
                                    length_afs[i].append(l)
                                    break
    
                    print(sum([len(b) for b in length_afs]))
                    # create boxplots
                    fig = plt.figure()
                    ax = plt.axes()		
                    bp = ax.boxplot(length_afs)
                    length_labels = []
                    for i,x in enumerate(length_intervals):
                        label = '<=' + str(i+1) + '/' + str(n_intervals)
                        length_labels.append(label)
    #				ax.set_yticks([0] + [i[1] for i in length_intervals])
    #				ax.set_yticklabels(["0"] + [str(i) + '/' + str(n_intervals) for i in range(1,n_intervals+1)])
                    ax.set_xticklabels(length_labels, rotation=45)
            #		ax.set_yscale("log")
                    plt.tight_layout()
                    pdf.savefig()
                    plt.close()


    for variant in variant_types.keys():
        if 'large' in variant:
            column_name = 'score_' + variant
            df['score_SVR'].fillna(df[column_name], inplace=True)
    
    df['confidence_level'] = 0
    df.confidence_level.where( (df.score_SVR.isnull() | (df.score_SVR<-0.5)), 1, inplace=True )
    df.confidence_level.where( (df.score_SVR.isnull() | (df.score_SVR<0.0)), 2, inplace=True )
    df.confidence_level.where( (df.score_SVR.isnull() | (df.score_SVR<0.5)), 3, inplace=True )
    df.confidence_level.where(~df.all_pass, 4, inplace=True)

    df = df.assign(is_large_insertion = lambda df: df.VARIANT_ID.isin(large_insertions))
    df = df.assign(is_large_deletion = lambda df: df.VARIANT_ID.isin(large_deletions))
    df = df.assign(is_large_complex = lambda df: df.VARIANT_ID.isin(large_complex))

    header = ["VARIANT_ID", "ac0_fail", "mendel_fail", "maxr2_fail", "hwe_fail", "carrier_gq_fail", "non_carrier_gq_fail", "score_SVR", "all_pass", "confidence_level", "is_large_insertion", "is_large_deletion", "is_large_complex"]
    df.to_csv(outname + '_filters.tsv', sep='\t', index=False, na_rep='nan')
