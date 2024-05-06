import argparse
from cyvcf2 import VCF
import numpy as np
import copy
import sys
from scipy import stats
from collections import namedtuple
from collections import defaultdict

# No longer needed. Using the Exact Test value from bcftools CLI
def check_HWE_chi(aaf, hom_ref, het, hom_alt, tot):
    
    if aaf == 0:
        assert het == 0
        assert hom_alt == 0
        return ("NaN", "NaN")
    
    exp_het = 2*aaf*(1-aaf)*tot
    exp_hom_ref = ((1-aaf)**2)*tot
    exp_hom_alt = (aaf**2)*tot
    try:
        chi_het = ((het-exp_het)**2)/exp_het
    except ZeroDivisionError:
        chi_het = 1e-20
    try:
        chi_alt = ((hom_alt-exp_hom_alt)**2)/exp_hom_alt
    except ZeroDivisionError:
        chi_alt = 1e-20
    try:
        chi_ref = ((hom_ref-exp_hom_ref)**2)/exp_hom_ref
    except ZeroDivisionError:
        chi_ref = 1e-20
    
    chi = chi_het+chi_ref+chi_alt
    return (np.round(chi,10), np.round(1-stats.chi2.cdf(chi, 1), 10))

def compute_panel_allele_statistics(record):
    """
    Compute allele related statistics.
    """
    an = 0
    ac = 0
    ac_aug = 0
    unknown = 0
    for genotype in record.genotypes:
        alleles = genotype[:-1]
        assert 1 <= len(alleles) <= 2
        # pseudohaplotypes
        if len(alleles) == 1:
            if alleles[0] != -1 and alleles[0] != 0:
                ac_aug += 1
        for a in alleles:
            if a == -1:
                unknown += 1
                continue
            assert a in [0,1]
            an += 1
            ac += a
    if an < 1:
        assert ac < 1
    af = ac / max(1.0, float(an))

    return (ac, ac_aug, an, unknown, af)
    
def compute_callset_statistics(record):
    """
    Compute genotype related statistics.
    """
    
    het_genotypes = 0
    hom_ref_genotypes = 0
    hom_alt_genotypes = 0
    total_genotypes = 0
    
    an = 0
    ac = 0
    af = 0
    unknown = 0

    for genotype in record.genotypes:
        alleles = genotype[:-1]
        assert 1 <= len(alleles) <= 2
        
        # this is only possible for unknown genotypes
        if len(alleles) == 1:
            # haploid genotype
            alleles.append(alleles[0])
        
        for a in alleles:
            if a == -1:
                unknown += 1
                continue
            assert a in [0, 1]
            an += 1
            ac += a
            
        if not -1 in alleles:
            total_genotypes += 1
            if sum(alleles) == 0:
                assert alleles == [0,0]
                hom_ref_genotypes += 1
            elif sum(alleles) == 1:
                assert 0 in alleles
                assert 1 in alleles
                het_genotypes += 1
            elif sum(alleles) == 2:
                assert alleles == [1,1]
                hom_alt_genotypes += 1
            else:
                sys.stderr.write("Inconsistent allele detected for %s"%(record.INFO['ID']))
    
    af = ac/max(1.0, float(an))
    allele_stats = (ac, af, an, unknown)
    
   
    genotype_stats = (hom_ref_genotypes, het_genotypes, hom_alt_genotypes, total_genotypes)
    
    hwe_stats = (np.round(record.INFO['HWE'],10), *check_HWE_chi(af, hom_ref_genotypes, het_genotypes, hom_alt_genotypes, total_genotypes))

    return allele_stats, genotype_stats, hwe_stats

def read_bubble_data(panel):
    
    panel = VCF(panel)
    # mapping variant IDs to bubble IDs
    var_id_to_bub = {}
    bub_id_to_allele_len = {}
    for variant in panel:
        bub_id = variant.ID
        var_ids = variant.INFO['ID']
        ref_length = len(variant.REF)
        alt_lengths = ','.join([str(len(x)) for x in variant.ALT])
        bub_id_to_allele_len[bub_id] = [ref_length, alt_lengths]
        for var_id in var_ids.split(','):
            for v in var_id.split(':'):
                if v in var_id_to_bub:
                    assert bub_id == var_id_to_bub[v]
                var_id_to_bub[v] = bub_id
    
    return var_id_to_bub, bub_id_to_allele_len

def get_sample_index_for_families(families, samples):
    family_to_index = {}
    for family in families:
        family_to_index[family] = []
        for sample in families[family]:
            family_to_index[family].append(samples.index(sample))
    return family_to_index


parser = argparse.ArgumentParser(prog='collect-vcf-stats.py', description="Collects the stats from the giggles callset vcf")
parser.add_argument('-bi-panel', help='Biallelic panel VCF.')
parser.add_argument('-bi-callset', help='Giggles multisample biallelic VCF')
parser.add_argument('-multi-panel', help='Multiallelic panel VCF')
parser.add_argument('-sample-sheet', metavar='sample_sheet', help='sample sheet containing list of unrelated samples')
args = parser.parse_args()

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

# reading bubble data and mapping variant ids to its bubbles
var_id_to_bub, bub_id_to_allele_len = read_bubble_data(args.multi_panel)

# header for variant info and bubble infor
header = ['VARIANT_ID',
        'VARIANT_TYPE',
        'VARIANT_LENGTH',
        'BUBBLE_ID',
        'BUBBLE_NUM_PATHS',
        'BUBBLE_ALT_PATHS_LENGTH',
        'BUBBLE_REF_PATH_LENGTH']
# header for panel info
header.extend(['PANEL_AC',
        'PANEL_AF',
        'PANEL_AC_PSEUDO',
        'PANEL_AN',
        'PANEL_NUM_UNKNOWN_HAPLOTYPES'])

# header for unrelated sample callset info
header.extend(['CALLSET_UNRELATED_AC', 
        'CALLSET_UNRELATED_AF',
        'CALLSET_UNRELATED_AN',
        'CALLSET_UNRELATED_NUM_UNTYPED_GENOTYPES',
        'CALLSET_UNRELATED_NUM_HET_SAMPLES',
        'CALLSET_UNRELATED_NUM_HOM_ALT_SAMPLES',
        'CALLSET_UNRELATED_NUM_HOM_REF_SAMPLES',
        'CALLSET_UNRELATED_NUM_TOTAL_SAMPLES'])

# header for HWE stats
header.extend(['HWE_EXACTTEST_STATISTIC',
        'HWE_CHI_TEST_STAT',
        'HWE_CHI_P_VALUE'])


empty_line = {key: None for key in header}
table_data = {}

panel_reader = VCF(args.bi_panel)
panel_samples = list(panel_reader.samples)

callset_reader = VCF(args.bi_callset)
callset_samples = list(callset_reader.samples)

panel_HG01258_gts = {}

print('\t'.join(header))

for variant in panel_reader:
    # require bi-allelic vcf with IDs
    assert len(variant.ALT) == 1
    var_id = variant.INFO['ID']
    ac, ac_aug, an, unknown, af = compute_panel_allele_statistics(variant)
    line = copy.deepcopy(empty_line)

    _,_,var_type,_,var_len = var_id.split('-')
    line['VARIANT_ID'] = var_id
    line['VARIANT_TYPE'] = var_type
    line['VARIANT_LENGTH'] = var_len

    bub_id = var_id_to_bub[var_id]
    num_paths = len(bub_id_to_allele_len[bub_id][1].split(',')) + 1
    line['BUBBLE_ID'] = bub_id
    line['BUBBLE_NUM_PATHS'] = num_paths
    line['BUBBLE_ALT_PATHS_LENGTH'] = bub_id_to_allele_len[bub_id][1]
    line['BUBBLE_REF_PATH_LENGTH'] = bub_id_to_allele_len[bub_id][0]

    line['PANEL_AC'] = ac
    line['PANEL_AF'] = af
    line['PANEL_AC_PSEUDO'] = ac_aug
    line['PANEL_AN'] = an
    line['PANEL_NUM_UNKNOWN_HAPLOTYPES'] = unknown
    table_data[var_id] = line

for n, variant in enumerate(callset_reader):
    assert len(variant.ALT) == 1
    var_id = variant.ID
    unrelated_allele_stats, unrelated_genotype_stats, hwe_stats = compute_callset_statistics(variant)

    table_data[var_id]['CALLSET_UNRELATED_AC'] = unrelated_allele_stats[0]
    table_data[var_id]['CALLSET_UNRELATED_AF'] = unrelated_allele_stats[1]
    table_data[var_id]['CALLSET_UNRELATED_AN'] = unrelated_allele_stats[2]
    table_data[var_id]['CALLSET_UNRELATED_NUM_UNTYPED_GENOTYPES'] = unrelated_allele_stats[3]
    table_data[var_id]['CALLSET_UNRELATED_NUM_HET_SAMPLES'] = unrelated_genotype_stats[1]
    table_data[var_id]['CALLSET_UNRELATED_NUM_HOM_ALT_SAMPLES'] = unrelated_genotype_stats[2]
    table_data[var_id]['CALLSET_UNRELATED_NUM_HOM_REF_SAMPLES'] = unrelated_genotype_stats[0]
    table_data[var_id]['CALLSET_UNRELATED_NUM_TOTAL_SAMPLES'] = unrelated_genotype_stats[3]

    table_data[var_id]['HWE_EXACTTEST_STATISTIC'] = hwe_stats[0]
    table_data[var_id]['HWE_CHI_TEST_STAT'] = hwe_stats[1]
    table_data[var_id]['HWE_CHI_P_VALUE'] = hwe_stats[2]

    line = '\t'.join([str(table_data[var_id][key]) for key in header])
    print(line)
    if (n)%1000 == 0:
        sys.stderr.write("\tRead %d variant records.\n"%(n))

sys.stderr.write("Completed generating callset statistics.\n")

