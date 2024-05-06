'''
Make the supplementary table that states which sample belongs to which list
'''

vienna_file = '/home/samarendra/Work/genotyping-results-localcopy/1kg-vienna-sample-list.tsv'
nygc_file = '/home/samarendra/Work/genotyping-results-localcopy/3202-nygc-sample-list.tsv'
kgp_file = '/home/samarendra/Work/genotyping-results-localcopy/2504-unrelated-sample-list.tsv'
hprc_file = '/home/samarendra/Work/genotyping-results-localcopy/47-hprc-sample-list.tsv'
hgsvc_file = '/home/samarendra/Work/genotyping-results-localcopy/44-hgsvc-phase2-sample-list.tsv'
hgsvc3_file = '/home/samarendra/Work/genotyping-results-localcopy/68-hgsvc-phase3-sample-list.tsv'
no_new_variant_samples = ['HG02661', 'HG01600', 'HG01308', 'NA19676']

def read_file(file):
    reader = open(file, 'r')
    samples = []
    for line in reader:
        samples.append(line.rstrip())
    reader.close()
    
    return samples

vienna_samples = read_file(vienna_file)
nygc_samples = read_file(nygc_file)
kgp_samples = read_file(kgp_file)
hprc_samples = read_file(hprc_file)
hgsvc_samples = read_file(hgsvc_file)
hgsvc3_samples = read_file(hgsvc3_file)

header = [
    'Sample',
    'in_1000GP_ONT',
    'in_NYGC',
    'in_1000GP_Phase3',
    'in_HPRC_Yr1',
    'in_HGSVC_Phase2',
    'in_HGSVC_Phase3',
    'by_SVarp',
    'by_Sniffles',
    'by_CuteSV',
    'by_Delly',
    'by_Giggles'
]

print('\t'.join(header))

for sample in vienna_samples:
    line = ""
    line += sample+'\t1\t'
    in_nygc = None
    in_kgp = None
    in_hprc = None
    in_hgsvc = None
    in_hgsvc3 = None
    by_svarp = None
    by_sniffles = '1'
    by_cutesv = '1'
    by_delly = '1'
    by_giggles = None
    
    in_nygc = '1' if sample in nygc_samples else '0'
    in_kgp = '1' if sample in kgp_samples else '0'
    in_hprc = '1' if sample in hprc_samples else '0'
    in_hgsvc = '1' if sample in hgsvc_samples else '0'
    in_hgsvc3 = '1' if sample in hgsvc3_samples else '0'

    by_svarp = '1' if (in_nygc == '1' and sample not in no_new_variant_samples) else '0'
    by_giggles = '1' if in_nygc == '1' else '0'

    line += '\t'.join([
        in_nygc,
        in_kgp,
        in_hprc,
        in_hgsvc,
        in_hgsvc3,
        by_svarp,
        by_sniffles,
        by_cutesv,
        by_delly,
        by_giggles])
    
    print(line)

