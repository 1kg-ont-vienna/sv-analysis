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
    

def compute_callset_statistics(record, qualities):
    """
    Compute genotype related statistics.
    """
    
    counts = defaultdict(int)
    het_genotypes = 0
    hom_ref_genotypes = 0
    hom_alt_genotypes = 0
    total_genotypes = 0
    
    an = 0
    ac = 0
    af = 0
    unknown = 0
    total_gq = 0

    gqs = record.format('GQ')
    total_gq = sum([q[0] for q in gqs])
    carrier_gq = 0
    non_carrier_gq = 0
    if total_gq < 0: total_gq=0

    for genotype, quality in zip(record.genotypes, gqs):
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
                non_carrier_gq += quality[0]
            elif sum(alleles) == 1:
                assert 0 in alleles
                assert 1 in alleles
                het_genotypes += 1
                carrier_gq += quality[0]
            elif sum(alleles) == 2:
                assert alleles == [1,1]
                hom_alt_genotypes += 1
                carrier_gq += quality[0]
            else:
                sys.stderr.write("Inconsistent allele detected for %s"%(record.INFO['ID']))
            for q in qualities:
                if quality[0] >= q:
                    counts[q] += 1
    
    af = ac/max(1.0, float(an))
    allele_stats = (ac, af, an, unknown)
    
    total_gq = total_gq/total_genotypes if total_genotypes != 0 else 0
    carrier_gq = carrier_gq/(het_genotypes + hom_alt_genotypes) if ((het_genotypes + hom_alt_genotypes) != 0) else 0
    non_carrier_gq = non_carrier_gq/(hom_ref_genotypes) if (hom_ref_genotypes != 0) else 0
    genotype_stats = (hom_ref_genotypes, het_genotypes, hom_alt_genotypes, total_genotypes, counts, total_gq, carrier_gq, non_carrier_gq)
    
    hwe_stats = (np.round(record.INFO['HWE'],10), *check_HWE_chi(af, hom_ref_genotypes, het_genotypes, hom_alt_genotypes, total_genotypes))

    return allele_stats, genotype_stats, hwe_stats

def compute_unrelated_callset_statistics(record, qualities, samples, unrelated_map):
    """
    Compute genotype related statistics.
    """
    
    counts = defaultdict(int)
    het_genotypes = 0
    hom_ref_genotypes = 0
    hom_alt_genotypes = 0
    total_genotypes = 0
    
    an = 0
    ac = 0
    af = 0
    unknown = 0
    total_gq = 0

    gqs = record.format('GQ')
    total_gq = sum([q[0] for q in gqs])
    carrier_gq = 0
    non_carrier_gq = 0
    if total_gq < 0: total_gq=0

    for genotype, quality, sample in zip(record.genotypes, gqs, samples):
        if unrelated_map[sample] == False:
            continue
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
                non_carrier_gq += quality[0]
            elif sum(alleles) == 1:
                assert 0 in alleles
                assert 1 in alleles
                het_genotypes += 1
                carrier_gq += quality[0]
            elif sum(alleles) == 2:
                assert alleles == [1,1]
                hom_alt_genotypes += 1
                carrier_gq += quality[0]
            else:
                sys.stderr.write("Inconsistent allele detected for %s"%(record.INFO['ID']))
            for q in qualities:
                if quality[0] >= q:
                    counts[q] += 1
    
    af = ac/max(1.0, float(an))
    allele_stats = (ac, af, an, unknown)
    
    total_gq = total_gq/total_genotypes if total_genotypes != 0 else 0
    carrier_gq = carrier_gq/(het_genotypes + hom_alt_genotypes) if ((het_genotypes + hom_alt_genotypes) != 0) else 0
    non_carrier_gq = non_carrier_gq/(hom_ref_genotypes) if (hom_ref_genotypes != 0) else 0
    genotype_stats = (hom_ref_genotypes, het_genotypes, hom_alt_genotypes, total_genotypes, counts, total_gq, carrier_gq, non_carrier_gq)

    return allele_stats, genotype_stats

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

class Genotype:
    def __init__(self, alleles, is_phased):
        self._alleles = alleles
        self._phased = is_phased

    def __str__(self):
        if self._phased:
            return '|'.join([str(i) for i in self._alleles])
        else:
            return '/'.join([str(i) for i in self._alleles])

    def get_alleles(self):
        return self._alleles
    
    def get_index(self):
        return sum(self._alleles)

    def get_ploidy(self):
        return len(self._alleles)

    def is_phased(self):
        return self._phased

    def __eq__(self, other):
        return sorted(self._alleles) == sorted(other._alleles)

    def is_hom_ref(self):
        return all([int(a) == 0 for a in self._alleles])

    def is_none(self):
        return self._alleles == []

def genotype_from_string(gt):
    is_phased = gt[-1]
    alleles = []
    if -1 in gt:
        # untyped
        return Genotype(alleles, is_phased)
    alleles = [int(allele) for allele in gt[:-1]]
    return Genotype(alleles, is_phased)

def get_sample_index_for_families(families, samples):
    family_to_index = {}
    for family in families:
        family_to_index[family] = []
        for sample in families[family]:
            family_to_index[family].append(samples.index(sample))
    return family_to_index

def check_mendelian_consistency(child_gt, parent1_gt, parent2_gt):
    child = child_gt.get_alleles()
    parent1 = parent1_gt.get_alleles()
    parent2 = parent2_gt.get_alleles()
    if child[0] in parent1 and child[1] in parent2:
        return True
    if child[0] in parent2 and child[1] in parent1:
        return True
    return False

def get_mendelian_statistics(record, family_to_index):
    gts = record.genotypes
    total_consistent_trios = 0
    alt_transmitted = 0
    all_het = 0
    all_abs = 0
    all_present = 0
    consistent_trios = 0
    inconsistent_trios = 0
    untyped_trios = 0
    total_untrivial_trios = 0
    
    for family_name, indexes in family_to_index.items():

        paternal_index, maternal_index, child_index = indexes
        gt_child = genotype_from_string(gts[child_index])
        gt_father = genotype_from_string(gts[paternal_index])
        gt_mother = genotype_from_string(gts[maternal_index])

        if any([g.is_none() for g in [gt_child, gt_father, gt_mother]]):
            untyped_trios += 1
        elif all([ g == gt_child for g in [gt_father, gt_mother]]):
            # all genotypes same, automatically consistent
            if gt_child == Genotype([0,0], False):
                all_abs += 1
            elif gt_child == Genotype([0,1], False):
                all_het += 1
            elif gt_child == Genotype([1,1], False):
                all_present += 1
            else:
                assert(False)
            total_consistent_trios += 1
            
        else:
            total_untrivial_trios += 1
            consistent = check_mendelian_consistency(gt_child, gt_father, gt_mother)
            if consistent:
                consistent_trios += 1
                total_consistent_trios += 1
                # check how often alt allele was transmitted
                alt_transmitted += sum(a!=0 for a in gt_child.get_alleles())
            else:
                inconsistent_trios += 1

    return (all_abs, all_het, all_present, consistent_trios, inconsistent_trios, untyped_trios, total_untrivial_trios, alt_transmitted)

def extract_genotype(record, index):
    gt = record.genotypes[index][:-1]
    if any([g == -1 for g in gt]):
        return -1
    return sum(gt)

def compare_genotypes(callset_gt, panel_gt):
    p_gt = "NaN" if panel_gt == -1 else panel_gt
    c_gt = "NaN" if callset_gt == -1 else callset_gt
    if (p_gt == 'NaN' or c_gt == 'NaN'):
        is_conc = 'NaN'
    else:
        is_conc = 1 if (p_gt == c_gt) else 0    
    
    return p_gt, c_gt, is_conc


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

qualities = [20, 50, 100, 250]

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
# header for callset info
header.extend(['CALLSET_AC', 
        'CALLSET_AF',
        'CALLSET_AN',
        'CALLSET_NUM_UNTYPED_GENOTYPES',
        'CALLSET_NUM_HET_SAMPLES',
        'CALLSET_NUM_HOM_ALT_SAMPLES',
        'CALLSET_NUM_HOM_REF_SAMPLES',
        'CALLSET_NUM_TOTAL_SAMPLES'])
for q in qualities:
        header.append('CALLSET_NUM_SAMPLES_GQ>=%s'%(str(q)))
header.extend(['CALLSET_AVG_GQ', 'CALLSET_CARRIER_AVG_GQ', 'CALLSET_NON_CARRIER_AVG_GQ'])

# header for unrelated sample callset info
header.extend(['CALLSET_UNRELATED_AC', 
        'CALLSET_UNRELATED_AF',
        'CALLSET_UNRELATED_AN',
        'CALLSET_UNRELATED_NUM_UNTYPED_GENOTYPES',
        'CALLSET_UNRELATED_NUM_HET_SAMPLES',
        'CALLSET_UNRELATED_NUM_HOM_ALT_SAMPLES',
        'CALLSET_UNRELATED_NUM_HOM_REF_SAMPLES',
        'CALLSET_UNRELATED_NUM_TOTAL_SAMPLES'])
for q in qualities:
        header.append('CALLSET_UNRELATED_NUM_SAMPLES_GQ>=%s'%(str(q)))
header.extend(['CALLSET_UNRELATED_AVG_GQ', 'CALLSET_UNRELATED_CARRIER_AVG_GQ', 'CALLSET_UNRELATED_NON_CARRIER_AVG_GQ'])

# header for HWE stats
header.extend(['HWE_EXACTTEST_STATISTIC',
        'HWE_CHI_TEST_STAT',
        'HWE_CHI_P_VALUE'])

# header for mendelian stats
header.extend(['MENDEL_ALL_ABSENT',
        'MENDEL_ALL_HET',
        'MENDEL_ALL_PRESENT',
        'MENDEL_UNTRIVIAL_CONSISTENT_TRIO',
        'MENDEL_INCONSISTENT_TRIO',
        'MENDEL_UNTYPED_TRIO',
        'MENDEL_TOTAL_UNTRIVIAL_TRIO',
        'MENDEL_ALT_TRANSMITTED'])
# header for self genotype concordance
header.extend(['CONCORDANCE_PANEL_GT',
        'CONCORDANCE_CALLSET_GT',
        'CONCORDANCE_IS_CONC'])


empty_line = {key: None for key in header}
table_data = {}

panel_reader = VCF(args.bi_panel)
panel_samples = list(panel_reader.samples)
panel_HG01258_index = panel_samples.index('HG01258')

callset_reader = VCF(args.bi_callset)
callset_samples = list(callset_reader.samples)
callset_HG01258_index = callset_samples.index('HG01258')

families = {
  '2418': 'NA19818_NA19819_NA19828'.split('_'),
  'CLM16': 'HG01256_HG01257_HG01258'.split('_'),
  'SH006': 'HG00418_HG00419_HG00420'.split('_'),
  'Y077': 'NA19128_NA19127_NA19129'.split('_'),
  '1463_Paternal': 'NA12889_NA12890_NA12877'.split('_'),
  '1463_Maternal': 'NA12891_NA12892_NA12878'.split('_')
}

family_to_index = get_sample_index_for_families(families, list(callset_reader.samples))

panel_HG01258_gts = {}

print('\t'.join(header))

for variant in panel_reader:
    # require bi-allelic vcf with IDs
    assert len(variant.ALT) == 1
    var_id = variant.INFO['ID']
    ac, ac_aug, an, unknown, af = compute_panel_allele_statistics(variant)
    line = copy.deepcopy(empty_line)
    
    panel_HG01258_gts[var_id] = extract_genotype(variant, panel_HG01258_index)
    
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
    var_id = variant.INFO['ID']
    allele_stats, genotype_stats, hwe_stats = compute_callset_statistics(variant, qualities)
    unrelated_allele_stats, unrelated_genotype_stats = compute_unrelated_callset_statistics(variant, qualities, callset_samples, unrelated_map)
    if var_id.split('-')[0] == 'chrX':
        mendel_stats = ["NaN"]*8
    else:
        mendel_stats = get_mendelian_statistics(variant, family_to_index)
    
    panel_gt, callset_gt, is_concordant = compare_genotypes(extract_genotype(variant, callset_HG01258_index), panel_HG01258_gts[var_id])

    table_data[var_id]['CALLSET_AC'] = allele_stats[0]
    table_data[var_id]['CALLSET_AF'] = allele_stats[1]
    table_data[var_id]['CALLSET_AN'] = allele_stats[2]
    table_data[var_id]['CALLSET_NUM_UNTYPED_GENOTYPES'] = allele_stats[3]
    table_data[var_id]['CALLSET_NUM_HET_SAMPLES'] = genotype_stats[1]
    table_data[var_id]['CALLSET_NUM_HOM_ALT_SAMPLES'] = genotype_stats[2]
    table_data[var_id]['CALLSET_NUM_HOM_REF_SAMPLES'] = genotype_stats[0]
    table_data[var_id]['CALLSET_NUM_TOTAL_SAMPLES'] = genotype_stats[3]
    for q in qualities:
        table_data[var_id]['CALLSET_NUM_SAMPLES_GQ>=%s'%(str(q))] = genotype_stats[4][q]
    table_data[var_id]['CALLSET_AVG_GQ'] = genotype_stats[5]
    table_data[var_id]['CALLSET_CARRIER_AVG_GQ'] = genotype_stats[6]
    table_data[var_id]['CALLSET_NON_CARRIER_AVG_GQ'] = genotype_stats[7]

    table_data[var_id]['CALLSET_UNRELATED_AC'] = unrelated_allele_stats[0]
    table_data[var_id]['CALLSET_UNRELATED_AF'] = unrelated_allele_stats[1]
    table_data[var_id]['CALLSET_UNRELATED_AN'] = unrelated_allele_stats[2]
    table_data[var_id]['CALLSET_UNRELATED_NUM_UNTYPED_GENOTYPES'] = unrelated_allele_stats[3]
    table_data[var_id]['CALLSET_UNRELATED_NUM_HET_SAMPLES'] = unrelated_genotype_stats[1]
    table_data[var_id]['CALLSET_UNRELATED_NUM_HOM_ALT_SAMPLES'] = unrelated_genotype_stats[2]
    table_data[var_id]['CALLSET_UNRELATED_NUM_HOM_REF_SAMPLES'] = unrelated_genotype_stats[0]
    table_data[var_id]['CALLSET_UNRELATED_NUM_TOTAL_SAMPLES'] = unrelated_genotype_stats[3]
    for q in qualities:
        table_data[var_id]['CALLSET_UNRELATED_NUM_SAMPLES_GQ>=%s'%(str(q))] = unrelated_genotype_stats[4][q]
    table_data[var_id]['CALLSET_UNRELATED_AVG_GQ'] = unrelated_genotype_stats[5]
    table_data[var_id]['CALLSET_UNRELATED_CARRIER_AVG_GQ'] = unrelated_genotype_stats[6]
    table_data[var_id]['CALLSET_UNRELATED_NON_CARRIER_AVG_GQ'] = unrelated_genotype_stats[7]

    table_data[var_id]['HWE_EXACTTEST_STATISTIC'] = hwe_stats[0]
    table_data[var_id]['HWE_CHI_TEST_STAT'] = hwe_stats[1]
    table_data[var_id]['HWE_CHI_P_VALUE'] = hwe_stats[2]

    table_data[var_id]['MENDEL_ALL_ABSENT'] = mendel_stats[0]
    table_data[var_id]['MENDEL_ALL_HET'] = mendel_stats[1]
    table_data[var_id]['MENDEL_ALL_PRESENT'] = mendel_stats[2]
    table_data[var_id]['MENDEL_UNTRIVIAL_CONSISTENT_TRIO'] = mendel_stats[3]
    table_data[var_id]['MENDEL_INCONSISTENT_TRIO'] = mendel_stats[4]
    table_data[var_id]['MENDEL_UNTYPED_TRIO'] = mendel_stats[5]
    table_data[var_id]['MENDEL_TOTAL_UNTRIVIAL_TRIO'] = mendel_stats[6]
    table_data[var_id]['MENDEL_ALT_TRANSMITTED'] = mendel_stats[7]
    
    table_data[var_id]['CONCORDANCE_PANEL_GT'] = panel_gt
    table_data[var_id]['CONCORDANCE_CALLSET_GT'] = callset_gt
    table_data[var_id]['CONCORDANCE_IS_CONC'] = is_concordant

    line = '\t'.join([str(table_data[var_id][key]) for key in header])
    print(line)
    if (n)%1000 == 0:
        sys.stderr.write("\tRead %d variant records.\n"%(n))

sys.stderr.write("Completed generating callset statistics.\n")

