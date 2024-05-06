import sys
import pandas as pd
import gzip

if __name__ == "__main__":
    filename = sys.argv[1]
    inp_vcf = sys.argv[2]
    outname = sys.argv[3]

    # creating the final table with the filters

    df = pd.read_csv(filename, sep='\t', header=0)
    df = df.assign(mendelian_consistency=lambda d: d['MENDEL_INCONSISTENT_TRIO'] / d['MENDEL_TOTAL_UNTRIVIAL_TRIO'])
    df = df.assign(CALLSET_HETEROZYGOSITY=lambda d: d['CALLSET_NUM_HET_SAMPLES'] / d['CALLSET_NUM_TOTAL_SAMPLES'])
    df = df.assign(CALLSET_UNRELATED_HETEROZYGOSITY=lambda d: d['CALLSET_UNRELATED_NUM_HET_SAMPLES'] / d['CALLSET_UNRELATED_NUM_TOTAL_SAMPLES'])

    # consider allele frequency computed across all genotyped samples
    df = df.assign(ac0_fail = lambda df: ~(df['CALLSET_AC'] > 0))
    # fail if mendelian inconsistent for 2 families
    df = df.assign(mendel_fail = lambda df: ((df.mendelian_consistency > 0.1)))
    # fail if average gq < 50
    df = df.assign(gq_fail = lambda df: df['CALLSET_AVG_GQ'] < 200)
    # fail if variant was incorrectly self-genotyped in HG01258
    df = df.assign(self_fail = lambda df: (df['CONCORDANCE_IS_CONC'] != 1))
    # defining strict set
    df = df.assign(all_pass = lambda df: ~(df.ac0_fail | df.mendel_fail | df.gq_fail | df.self_fail))
    
    # new strict set
    #df = df.assign(new_pass = lambda df: (df.all_pass) | (df['HWE_EXACTTEST_STATISTIC'] > 0.5) | (df['CALLSET_AF'] > 0.5))
    #df = df.assign(new_pass = lambda df: (df.all_pass))

    df.to_csv(outname+'qc.table.final.tsv', sep='\t', index=False)

    passed_var_ids = df[df['new_pass'] == True]['VARIANT_ID'].to_numpy()


    # subsetting vcf
    out_vcf = open(outname+'filtered-final.vcf', 'w')
    out_vcf_lowq = open(outname+'lowq-final.vcf', 'w')
    vcf_reader = gzip.open(inp_vcf, 'r')
    count = 0
    for line in vcf_reader:
        line = line.decode('utf-8')
        if line[0] == '#':
            out_vcf.write(line)
            out_vcf_lowq.write(line)
            continue
        id = line.split('\t')[2]
        if count == passed_var_ids.__len__():
            break
        if id == passed_var_ids[count]:
            count += 1
            out_vcf.write(line)
        else:
            out_vcf_lowq.write(line)
