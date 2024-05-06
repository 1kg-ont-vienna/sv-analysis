import argparse
from cyvcf2 import VCF
import pandas
import sys
from collections import namedtuple
import copy

AlleleStats = namedtuple('AlleleStats','af ac an untyped')
GenotypeStats = namedtuple('GenotypeStats', 'heterozygosity het hom_ref hom_alt total chi passed')

def check_HWE_chi(aaf, hom_ref, het, hom_alt, tot):
    
    if aaf == 0:
        assert het == 0
        assert hom_alt == 0
        return ("NaN", 0)
    
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
    if chi < 3.841:
        return (chi, 1)
    else:
        return (chi, 0)

def create_GenotypeStats(ref, het, alt, tot, chi, passed):
    return GenotypeStats( str(het / max(1.0, float(tot))), str(het), str(ref), str(alt), str(tot), str(chi), str(passed))

def compute_panel_allele_statistics(record):
    """
    Compute allele related statistics.
    """
    an = 0
    ac = 0
    unknown = 0
    for genotype in record.genotypes:
        alleles = genotype[:-1]
        assert 1 <= len(alleles) <= 2
        
        #### TO KEEP OR NOT TO KEEP
        if len(alleles) == 1:
            # haploid genotype
            alleles.append(alleles[0])
        ####
        
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
    return AlleleStats(str(af), str(ac), str(an), str(unknown))

def compute_callset_statistics(record, sample_pop_map, samples, sample2coverage):
    """
    Compute genotype related statistics.
    """
    het_genotypes = {
        'high': {'all': 0, 'AFR': 0, 'AMR': 0, 'EAS': 0, 'EUR': 0, 'SAS': 0}, 
        'low': {'all': 0, 'AFR': 0, 'AMR': 0, 'EAS': 0, 'EUR': 0, 'SAS': 0},
        'all': {'all': 0, 'AFR': 0, 'AMR': 0, 'EAS': 0, 'EUR': 0, 'SAS': 0}
        }
    hom_ref_genotypes = copy.deepcopy(het_genotypes)
    hom_alt_genotypes = copy.deepcopy(het_genotypes)
    total_genotypes = copy.deepcopy(het_genotypes)
    
    an = copy.deepcopy(het_genotypes)
    ac = copy.deepcopy(het_genotypes)
    af = copy.deepcopy(het_genotypes)
    unknown = copy.deepcopy(het_genotypes)

    for genotype, sample in zip(record.genotypes, samples):
        pop_code = sample_pop_map[sample]
        cov = 'high' if sample2coverage[sample] else 'low'
        assert pop_code in ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']
        alleles = genotype[:-1]
        assert 1 <= len(alleles) <= 2
        if -1 in alleles:
            unknown[cov]['all'] += 2
            unknown[cov][pop_code] += 2
            unknown['all']['all'] += 2
            unknown['all'][pop_code] += 2
            continue
        
        for a in alleles:
            assert a in [0, 1]
            an[cov]['all'] += 1
            an[cov][pop_code] += 1
            ac[cov]['all'] += a
            ac[cov][pop_code] += a
            an['all']['all'] += 1
            an['all'][pop_code] += 1
            ac['all']['all'] += a
            ac['all'][pop_code] += a

        total_genotypes[cov]['all'] += 1
        total_genotypes[cov][pop_code] += 1
        total_genotypes['all']['all'] += 1
        total_genotypes['all'][pop_code] += 1
        if sum(alleles) == 0:
            assert alleles == [0,0]
            hom_ref_genotypes[cov]['all'] += 1
            hom_ref_genotypes[cov][pop_code] += 1
            hom_ref_genotypes['all']['all'] += 1
            hom_ref_genotypes['all'][pop_code] += 1
            continue
        if sum(alleles) == 1:
            assert 0 in alleles
            assert 1 in alleles
            het_genotypes[cov]['all'] += 1
            het_genotypes[cov][pop_code] += 1
            het_genotypes['all']['all'] += 1
            het_genotypes['all'][pop_code] += 1
            continue
        if sum(alleles) == 2:
            assert alleles == [1,1]
            hom_alt_genotypes[cov]['all'] += 1
            hom_alt_genotypes[cov][pop_code] += 1
            hom_alt_genotypes['all']['all'] += 1
            hom_alt_genotypes['all'][pop_code] += 1
            continue
        sys.stderr.write("Inconsistent allele detected for %s"%(record.INFO['ID']))
            
    low_cov_allele_stats = {}
    high_cov_allele_stats = {}
    all_cov_allele_stats = {}
    for i in an['high'].keys():
        if an['high'][i] < 1:
            assert ac['high'][i] < 1
        af['high'][i] = ac['high'][i] / max(1.0, float(an['high'][i]))
        high_cov_allele_stats[i] = AlleleStats(str(af['high'][i]), str(ac['high'][i]), str(an['high'][i]), str(unknown['high'][i]))
    for i in an['low'].keys():
        if an['low'][i] < 1:
            assert ac['low'][i] < 1
        af['low'][i] = ac['low'][i] / max(1.0, float(an['low'][i]))
        low_cov_allele_stats[i] = AlleleStats(str(af['low'][i]), str(ac['low'][i]), str(an['low'][i]), str(unknown['low'][i]))
    for i in an['all'].keys():
        if an['all'][i] < 1:
            assert ac['all'][i] < 1
        af['all'][i] = ac['all'][i] / max(1.0, float(an['all'][i]))
        all_cov_allele_stats[i] = AlleleStats(str(af['all'][i]), str(ac['all'][i]), str(an['all'][i]), str(unknown['all'][i]))
    
    low_cov_genotype_stats = {}
    high_cov_genotype_stats = {}
    all_cov_genotype_stats = {}
    for i in het_genotypes['high'].keys():
        high_cov_genotype_stats[i] = create_GenotypeStats(hom_ref_genotypes['high'][i], het_genotypes['high'][i], hom_alt_genotypes['high'][i], total_genotypes['high'][i], *check_HWE_chi(af['high'][i], hom_ref_genotypes['high'][i], het_genotypes['high'][i], hom_alt_genotypes['high'][i], total_genotypes['high'][i]))
    for i in het_genotypes['low'].keys():
        low_cov_genotype_stats[i] = create_GenotypeStats(hom_ref_genotypes['low'][i], het_genotypes['low'][i], hom_alt_genotypes['low'][i], total_genotypes['low'][i], *check_HWE_chi(af['low'][i], hom_ref_genotypes['low'][i], het_genotypes['low'][i], hom_alt_genotypes['low'][i], total_genotypes['low'][i]))
    for i in het_genotypes['all'].keys():
        all_cov_genotype_stats[i] = create_GenotypeStats(hom_ref_genotypes['all'][i], het_genotypes['all'][i], hom_alt_genotypes['all'][i], total_genotypes['all'][i], *check_HWE_chi(af['all'][i], hom_ref_genotypes['all'][i], het_genotypes['all'][i], hom_alt_genotypes['all'][i], total_genotypes['all'][i]))
    
    return low_cov_allele_stats, high_cov_allele_stats, all_cov_allele_stats, low_cov_genotype_stats, high_cov_genotype_stats, all_cov_genotype_stats

def compute_unrelated_callset_statistics(record, sample_pop_map, samples, sample2coverage, unrelated_map):
    """
    Compute genotype related statistics.
    """
    het_genotypes = {
        'high': {'all': 0, 'AFR': 0, 'AMR': 0, 'EAS': 0, 'EUR': 0, 'SAS': 0}, 
        'low': {'all': 0, 'AFR': 0, 'AMR': 0, 'EAS': 0, 'EUR': 0, 'SAS': 0},
        'all': {'all': 0, 'AFR': 0, 'AMR': 0, 'EAS': 0, 'EUR': 0, 'SAS': 0}
        }
    hom_ref_genotypes = copy.deepcopy(het_genotypes)
    hom_alt_genotypes = copy.deepcopy(het_genotypes)
    total_genotypes = copy.deepcopy(het_genotypes)
    
    an = copy.deepcopy(het_genotypes)
    ac = copy.deepcopy(het_genotypes)
    af = copy.deepcopy(het_genotypes)
    unknown = copy.deepcopy(het_genotypes)

    for genotype, sample in zip(record.genotypes, samples):
        if unrelated_map[sample] == False:
            continue
        pop_code = sample_pop_map[sample]
        cov = 'high' if sample2coverage[sample] else 'low'
        assert pop_code in ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']
        alleles = genotype[:-1]
        assert 1 <= len(alleles) <= 2
        if -1 in alleles:
            unknown[cov]['all'] += 2
            unknown[cov][pop_code] += 2
            unknown['all']['all'] += 2
            unknown['all'][pop_code] += 2
            continue
        
        for a in alleles:
            assert a in [0, 1]
            an[cov]['all'] += 1
            an[cov][pop_code] += 1
            ac[cov]['all'] += a
            ac[cov][pop_code] += a
            an['all']['all'] += 1
            an['all'][pop_code] += 1
            ac['all']['all'] += a
            ac['all'][pop_code] += a

        total_genotypes[cov]['all'] += 1
        total_genotypes[cov][pop_code] += 1
        total_genotypes['all']['all'] += 1
        total_genotypes['all'][pop_code] += 1
        if sum(alleles) == 0:
            assert alleles == [0,0]
            hom_ref_genotypes[cov]['all'] += 1
            hom_ref_genotypes[cov][pop_code] += 1
            hom_ref_genotypes['all']['all'] += 1
            hom_ref_genotypes['all'][pop_code] += 1
            continue
        if sum(alleles) == 1:
            assert 0 in alleles
            assert 1 in alleles
            het_genotypes[cov]['all'] += 1
            het_genotypes[cov][pop_code] += 1
            het_genotypes['all']['all'] += 1
            het_genotypes['all'][pop_code] += 1
            continue
        if sum(alleles) == 2:
            assert alleles == [1,1]
            hom_alt_genotypes[cov]['all'] += 1
            hom_alt_genotypes[cov][pop_code] += 1
            hom_alt_genotypes['all']['all'] += 1
            hom_alt_genotypes['all'][pop_code] += 1
            continue
        sys.stderr.write("Inconsistent allele detected for %s"%(record.INFO['ID']))
            
    low_cov_allele_stats = {}
    high_cov_allele_stats = {}
    all_cov_allele_stats = {}
    for i in an['high'].keys():
        if an['high'][i] < 1:
            assert ac['high'][i] < 1
        af['high'][i] = ac['high'][i] / max(1.0, float(an['high'][i]))
        high_cov_allele_stats[i] = AlleleStats(str(af['high'][i]), str(ac['high'][i]), str(an['high'][i]), str(unknown['high'][i]))
    for i in an['low'].keys():
        if an['low'][i] < 1:
            assert ac['low'][i] < 1
        af['low'][i] = ac['low'][i] / max(1.0, float(an['low'][i]))
        low_cov_allele_stats[i] = AlleleStats(str(af['low'][i]), str(ac['low'][i]), str(an['low'][i]), str(unknown['low'][i]))
    for i in an['all'].keys():
        if an['all'][i] < 1:
            assert ac['all'][i] < 1
        af['all'][i] = ac['all'][i] / max(1.0, float(an['all'][i]))
        all_cov_allele_stats[i] = AlleleStats(str(af['all'][i]), str(ac['all'][i]), str(an['all'][i]), str(unknown['all'][i]))
    
    low_cov_genotype_stats = {}
    high_cov_genotype_stats = {}
    all_cov_genotype_stats = {}
    for i in het_genotypes['high'].keys():
        high_cov_genotype_stats[i] = create_GenotypeStats(hom_ref_genotypes['high'][i], het_genotypes['high'][i], hom_alt_genotypes['high'][i], total_genotypes['high'][i], *check_HWE_chi(af['high'][i], hom_ref_genotypes['high'][i], het_genotypes['high'][i], hom_alt_genotypes['high'][i], total_genotypes['high'][i]))
    for i in het_genotypes['low'].keys():
        low_cov_genotype_stats[i] = create_GenotypeStats(hom_ref_genotypes['low'][i], het_genotypes['low'][i], hom_alt_genotypes['low'][i], total_genotypes['low'][i], *check_HWE_chi(af['low'][i], hom_ref_genotypes['low'][i], het_genotypes['low'][i], hom_alt_genotypes['low'][i], total_genotypes['low'][i]))
    for i in het_genotypes['all'].keys():
        all_cov_genotype_stats[i] = create_GenotypeStats(hom_ref_genotypes['all'][i], het_genotypes['all'][i], hom_alt_genotypes['all'][i], total_genotypes['all'][i], *check_HWE_chi(af['all'][i], hom_ref_genotypes['all'][i], het_genotypes['all'][i], hom_alt_genotypes['all'][i], total_genotypes['all'][i]))
    
    return low_cov_allele_stats, high_cov_allele_stats, all_cov_allele_stats, low_cov_genotype_stats, high_cov_genotype_stats, all_cov_genotype_stats

parser = argparse.ArgumentParser(prog='hwe-af-stats.py', description="Collects the stats from the giggles callset vcf population and variant type stratified")
parser.add_argument('-meta', metavar='meta', help='Sample population data')
parser.add_argument('-sample-sheet', metavar='sample_sheet', help='sample sheet containing list of unrelated samples')
parser.add_argument('-bi-panel', metavar='bi_panel', help='Biallelic panel VCF.')
parser.add_argument('-multi-panel', metavar='multi_panel', help='Multiallelic panel VCF.')
parser.add_argument('-callset', metavar='callset', help='Giggles multisample biallelic VCF')
parser.add_argument('-output', metavar='output', help='Output prefix for tsv files')
parser.add_argument('-coverages', metavar='coverages', help='TSV file with coverages')
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

coverages = pandas.read_csv(args.coverages, sep='\t', header=0)
coverages = coverages[["SAMPLE", "T2T_median_cov"]]
sample2coverage = {}
for line in coverages.to_numpy():
    if line[1] > 15:
        sample2coverage[line[0]] = True
    else:
        sample2coverage[line[0]] = False

metadata = pandas.read_csv(args.meta, sep='\t', header=0)
metadata = metadata[["Sample name", "Population code", "Superpopulation code"]]
sample2pop = {}
for line in metadata.to_numpy():
    sample2pop[line[0]] = line[2]

panel_reader = VCF(args.bi_panel)
panel_stats = {}
callset_reader = VCF(args.callset)
callset_samples = callset_reader.samples
low_cov_callset_stats = {}
high_cov_callset_stats = {}
all_cov_callset_stats = {}
unrelated_low_cov_callset_stats = {}
unrelated_high_cov_callset_stats = {}
unrelated_all_cov_callset_stats = {}

for variant in panel_reader:
    # require bi-allelic vcf with IDs
    assert len(variant.ALT) == 1
    var_id = variant.INFO['ID']
    allele_stats = compute_panel_allele_statistics(variant)
    panel_stats[var_id] = allele_stats

panel_reader.close()

sys.stderr.write("Completed generating panel statistics.\n")
sys.stderr.write("Found %d variant IDs.\n"%(len(panel_stats)))
sys.stderr.write("\nReading Callset VCF.\n")

for n, variant in enumerate(callset_reader):
    assert len(variant.ALT) == 1
    var_id = variant.INFO['ID']
    
    low_cov_allele_stats, high_cov_allele_stats, all_cov_allele_stats, low_cov_genotype_stats, high_cov_genotype_stats, all_cov_genotype_stats = compute_callset_statistics(variant, sample2pop, callset_samples, sample2coverage)
    low_cov_callset_stats[var_id] = [low_cov_allele_stats, low_cov_genotype_stats]
    high_cov_callset_stats[var_id] = [high_cov_allele_stats, high_cov_genotype_stats]
    all_cov_callset_stats[var_id] = [all_cov_allele_stats, all_cov_genotype_stats]
    
    low_cov_allele_stats, high_cov_allele_stats, all_cov_allele_stats, low_cov_genotype_stats, high_cov_genotype_stats, all_cov_genotype_stats = compute_unrelated_callset_statistics(variant, sample2pop, callset_samples, sample2coverage, unrelated_map)
    unrelated_low_cov_callset_stats[var_id] = [low_cov_allele_stats, low_cov_genotype_stats]
    unrelated_high_cov_callset_stats[var_id] = [high_cov_allele_stats, high_cov_genotype_stats]
    unrelated_all_cov_callset_stats[var_id] = [all_cov_allele_stats, all_cov_genotype_stats]
    if (n+1)%1000 == 0:
        sys.stderr.write("\tRead %d variant records.\n"%(n))

panel = VCF(args.multi_panel)
# mapping variant IDs to bubble IDs
var_id_to_bub = {}
bub_id_to_allele_len = {}
for variant in panel:
    bub_id = variant.ID
    var_ids = variant.INFO['ID']
    ref_length = len(variant.REF)
    alt_lengths = ','.join([str(len(x)) for x in variant.ALT])
    bub_id_to_allele_len[bub_id] = [str(ref_length), str(alt_lengths)]
    for var_id in var_ids.split(','):
        for v in var_id.split(':'):
            if v in var_id_to_bub:
                assert bub_id == var_id_to_bub[v]
            var_id_to_bub[v] = str(bub_id)
panel.close()

high_cov_writer = open(args.output+'-high_cov.tsv', 'w')
low_cov_writer = open(args.output+'-low_cov.tsv', 'w')
all_cov_writer = open(args.output+'-all_cov.tsv', 'w')
# print stats for all IDs in genotypes VCFin
header = ['variant_id',
        'variant_type',
        'variant_length',
        'bub_id',
        'ref_allele_len',
        'alt_allele_len',
        'panel_allele_freq',
        'panel_alternative_alleles',
        'panel_total_alleles',
        'panel_unknown_alleles']

for c in ['allpop', 'AFR', 'AMR', 'EAS', 'EUR', 'SAS']:
    for u in ['all', 'unrelated']:
        header.append('%s-%s_allele_freq'%(c,u))
        header.append('%s-%s_alternative_alleles'%(c,u))
        header.append('%s-%s_total_alleles'%(c,u))
        header.append('%s-%s_unknown_alleles'%(c,u))
        header.append('%s-%s_heterozygosity'%(c,u))
        header.append('%s-%s_heterozygous_genotypes'%(c,u))
        header.append('%s-%s_homozygous_reference_genotypes'%(c,u))
        header.append('%s-%s_homozygous_alternate_genotypes'%(c,u))
        header.append('%s-%s_total_genotypes'%(c,u))
        header.append('%s-%s_chi_stat'%(c,u))
        header.append('%s-%s_chi_passed'%(c,u))

print('\t'.join(header), file=high_cov_writer)
print('\t'.join(header), file=low_cov_writer)
print('\t'.join(header), file=all_cov_writer)

for var_id in high_cov_callset_stats:
    if not var_id in panel_stats:
        continue
    
    _,_,var_type,_,var_len = var_id.split('-')
    line = [var_id,
            var_type,
            var_len,
            var_id_to_bub[var_id],
            bub_id_to_allele_len[var_id_to_bub[var_id]][0],
            bub_id_to_allele_len[var_id_to_bub[var_id]][1],
            panel_stats[var_id].af,
            panel_stats[var_id].ac,
            panel_stats[var_id].an,
            panel_stats[var_id].untyped]
    for c in ['all', 'AFR', 'AMR', 'EAS', 'EUR', 'SAS']:
        line.extend(
            [high_cov_callset_stats[var_id][0][c].af,
            high_cov_callset_stats[var_id][0][c].ac,
            high_cov_callset_stats[var_id][0][c].an,
            high_cov_callset_stats[var_id][0][c].untyped,
            high_cov_callset_stats[var_id][1][c].heterozygosity,
            high_cov_callset_stats[var_id][1][c].het,
            high_cov_callset_stats[var_id][1][c].hom_ref,
            high_cov_callset_stats[var_id][1][c].hom_alt,
            high_cov_callset_stats[var_id][1][c].total,
            high_cov_callset_stats[var_id][1][c].chi,
            high_cov_callset_stats[var_id][1][c].passed
            ])
        line.extend(
            [unrelated_high_cov_callset_stats[var_id][0][c].af,
            unrelated_high_cov_callset_stats[var_id][0][c].ac,
            unrelated_high_cov_callset_stats[var_id][0][c].an,
            unrelated_high_cov_callset_stats[var_id][0][c].untyped,
            unrelated_high_cov_callset_stats[var_id][1][c].heterozygosity,
            unrelated_high_cov_callset_stats[var_id][1][c].het,
            unrelated_high_cov_callset_stats[var_id][1][c].hom_ref,
            unrelated_high_cov_callset_stats[var_id][1][c].hom_alt,
            unrelated_high_cov_callset_stats[var_id][1][c].total,
            unrelated_high_cov_callset_stats[var_id][1][c].chi,
            unrelated_high_cov_callset_stats[var_id][1][c].passed
            ])
    assert len(line) == len(header)
    print('\t'.join(line), file=high_cov_writer)

for var_id in low_cov_callset_stats:
    if not var_id in panel_stats:
        continue
    
    _,_,var_type,_,var_len = var_id.split('-')
    line = [var_id,
            var_type,
            var_len,
            var_id_to_bub[var_id],
            bub_id_to_allele_len[var_id_to_bub[var_id]][0],
            bub_id_to_allele_len[var_id_to_bub[var_id]][1],
            panel_stats[var_id].af,
            panel_stats[var_id].ac,
            panel_stats[var_id].an,
            panel_stats[var_id].untyped]
    for c in ['all', 'AFR', 'AMR', 'EAS', 'EUR', 'SAS']:
        line.extend(
            [low_cov_callset_stats[var_id][0][c].af,
            low_cov_callset_stats[var_id][0][c].ac,
            low_cov_callset_stats[var_id][0][c].an,
            low_cov_callset_stats[var_id][0][c].untyped,
            low_cov_callset_stats[var_id][1][c].heterozygosity,
            low_cov_callset_stats[var_id][1][c].het,
            low_cov_callset_stats[var_id][1][c].hom_ref,
            low_cov_callset_stats[var_id][1][c].hom_alt,
            low_cov_callset_stats[var_id][1][c].total,
            low_cov_callset_stats[var_id][1][c].chi,
            low_cov_callset_stats[var_id][1][c].passed
            ])
        line.extend(
            [unrelated_low_cov_callset_stats[var_id][0][c].af,
            unrelated_low_cov_callset_stats[var_id][0][c].ac,
            unrelated_low_cov_callset_stats[var_id][0][c].an,
            unrelated_low_cov_callset_stats[var_id][0][c].untyped,
            unrelated_low_cov_callset_stats[var_id][1][c].heterozygosity,
            unrelated_low_cov_callset_stats[var_id][1][c].het,
            unrelated_low_cov_callset_stats[var_id][1][c].hom_ref,
            unrelated_low_cov_callset_stats[var_id][1][c].hom_alt,
            unrelated_low_cov_callset_stats[var_id][1][c].total,
            unrelated_low_cov_callset_stats[var_id][1][c].chi,
            unrelated_low_cov_callset_stats[var_id][1][c].passed
            ])
    assert len(line) == len(header)
    print('\t'.join(line), file=low_cov_writer)

for var_id in all_cov_callset_stats:
    if not var_id in panel_stats:
        continue
    
    _,_,var_type,_,var_len = var_id.split('-')
    line = [var_id,
            var_type,
            var_len,
            var_id_to_bub[var_id],
            bub_id_to_allele_len[var_id_to_bub[var_id]][0],
            bub_id_to_allele_len[var_id_to_bub[var_id]][1],
            panel_stats[var_id].af,
            panel_stats[var_id].ac,
            panel_stats[var_id].an,
            panel_stats[var_id].untyped]
    for c in ['all', 'AFR', 'AMR', 'EAS', 'EUR', 'SAS']:
        line.extend(
            [all_cov_callset_stats[var_id][0][c].af,
            all_cov_callset_stats[var_id][0][c].ac,
            all_cov_callset_stats[var_id][0][c].an,
            all_cov_callset_stats[var_id][0][c].untyped,
            all_cov_callset_stats[var_id][1][c].heterozygosity,
            all_cov_callset_stats[var_id][1][c].het,
            all_cov_callset_stats[var_id][1][c].hom_ref,
            all_cov_callset_stats[var_id][1][c].hom_alt,
            all_cov_callset_stats[var_id][1][c].total,
            all_cov_callset_stats[var_id][1][c].chi,
            all_cov_callset_stats[var_id][1][c].passed
            ])
        line.extend(
            [unrelated_all_cov_callset_stats[var_id][0][c].af,
            unrelated_all_cov_callset_stats[var_id][0][c].ac,
            unrelated_all_cov_callset_stats[var_id][0][c].an,
            unrelated_all_cov_callset_stats[var_id][0][c].untyped,
            unrelated_all_cov_callset_stats[var_id][1][c].heterozygosity,
            unrelated_all_cov_callset_stats[var_id][1][c].het,
            unrelated_all_cov_callset_stats[var_id][1][c].hom_ref,
            unrelated_all_cov_callset_stats[var_id][1][c].hom_alt,
            unrelated_all_cov_callset_stats[var_id][1][c].total,
            unrelated_all_cov_callset_stats[var_id][1][c].chi,
            unrelated_all_cov_callset_stats[var_id][1][c].passed
            ])
    assert len(line) == len(header)
    print('\t'.join(line), file=all_cov_writer)

high_cov_writer.close()
low_cov_writer.close()
all_cov_writer.close()