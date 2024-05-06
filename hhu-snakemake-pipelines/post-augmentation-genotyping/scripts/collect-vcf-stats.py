import argparse
from cyvcf2 import VCF
import pandas
import sys
from collections import namedtuple
from collections import defaultdict

AlleleStats = namedtuple('AlleleStats','af ac an untyped')
GenotypeStats = namedtuple('GenotypeStats', 'heterozygosity het hom_ref hom_alt total')

def create_GenotypeStats(ref, het, alt, tot):
    return GenotypeStats( str(het / max(1.0, float(tot))), str(het), str(ref), str(alt), str(tot))

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
        if len(alleles) == 1:
            # haploid genotype
            alleles.append(alleles[0])
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

def compute_callset_statistics(record, qualities=None, sample_pop_map=None, samples=None):
    """
    Compute genotype related statistics.
    """
    counts = {'all': defaultdict(int), 'AFR': defaultdict(int), 'AMR': defaultdict(int), 'EAS': defaultdict(int), 'EUR': defaultdict(int), 'SAS': defaultdict(int)}
    het_genotypes = {'all': 0, 'AFR': 0, 'AMR': 0, 'EAS': 0, 'EUR': 0, 'SAS': 0}
    hom_ref_genotypes = {'all': 0, 'AFR': 0, 'AMR': 0, 'EAS': 0, 'EUR': 0, 'SAS': 0}
    hom_alt_genotypes = {'all': 0, 'AFR': 0, 'AMR': 0, 'EAS': 0, 'EUR': 0, 'SAS': 0}
    total_genotypes = {'all': 0, 'AFR': 0, 'AMR': 0, 'EAS': 0, 'EUR': 0, 'SAS': 0}
    
    an = {'all': 0, 'AFR': 0, 'AMR': 0, 'EAS': 0, 'EUR': 0, 'SAS': 0}
    ac = {'all': 0, 'AFR': 0, 'AMR': 0, 'EAS': 0, 'EUR': 0, 'SAS': 0}
    af = {'all': 0, 'AFR': 0, 'AMR': 0, 'EAS': 0, 'EUR': 0, 'SAS': 0}
    unknown = {'all': 0, 'AFR': 0, 'AMR': 0, 'EAS': 0, 'EUR': 0, 'SAS': 0}

    gqs = record.format('GQ') if qualities is not None else [None]*len(record.genotypes)
    for genotype, quality, sample in zip(record.genotypes, gqs, samples):
        pop_code = sample_pop_map[sample]
        assert pop_code in ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']
        alleles = genotype[:-1]
        assert 1 <= len(alleles) <= 2
        if len(alleles) == 1:
            # haploid genotype
            alleles.append(alleles[0])
        
        for a in alleles:
            if a == -1:
                unknown['all'] += 1
                unknown[pop_code] += 1
                continue
            assert a in [0, 1]
            an['all'] += 1
            an[pop_code] += 1
            ac['all'] += a
            ac[pop_code] += a
        
        if not -1 in alleles:
            total_genotypes['all'] += 1
            total_genotypes[pop_code] += 1
            if sum(alleles) == 0:
                assert alleles == [0,0]
                hom_ref_genotypes['all'] += 1
                hom_ref_genotypes[pop_code] += 1
            elif sum(alleles) == 1:
                assert 0 in alleles
                assert 1 in alleles
                het_genotypes['all'] += 1
                het_genotypes[pop_code] += 1
            elif sum(alleles) == 2:
                assert alleles == [1,1]
                hom_alt_genotypes['all'] += 1
                hom_alt_genotypes[pop_code] += 1
            else:
                sys.stderr.write("Inconsistent allele detected for %s"%(record.INFO['ID']))
            # read GQ
            if qualities is not None:
                for q in qualities:
                    if int(quality[0]) >= q:
                        counts['all'][q] += 1
                        counts[pop_code][q] += 1
    
    allele_stats = {}
    for i in an.keys():
        if an[i] < 1:
            assert ac[i] < 1
        af[i] = ac[i] / max(1.0, float(an[i]))
        allele_stats[i] = AlleleStats(str(af[i]), str(ac[i]), str(an[i]), str(unknown[i]))

    genotype_stats = {}
    for i in het_genotypes.keys():
        genotype_stats[i] = create_GenotypeStats(hom_ref_genotypes[i], het_genotypes[i], hom_alt_genotypes[i], total_genotypes[i])
    
    return allele_stats, genotype_stats, counts

parser = argparse.ArgumentParser(prog='collect-vcf-stats.py', description="Collects the stats from the giggles callset vcf")
parser.add_argument('-meta', metavar='meta', help='Sample population data')
parser.add_argument('-panel', metavar='panel', help='Biallelic panel VCF.')
parser.add_argument('-callset', metavar='callset', help='Giggles multisample biallelic VCF')
args = parser.parse_args()

# reading metadata (sample-to-population map and color coding)
metadata = pandas.read_csv(args.meta, sep='\t', header=0)
metadata = metadata[["Sample name", "Population code", "Superpopulation code"]]
sample2pop = {}
for line in metadata.to_numpy():
    sample2pop[line[0]] = line[2]

panel_reader = VCF(args.panel)
panel_stats = {}
callset_reader = VCF(args.callset)
callset_samples = callset_reader.samples
callset_stats = {}

for variant in panel_reader:
    # require bi-allelic vcf with IDs
    assert len(variant.ALT) == 1
    var_id = variant.INFO['ID']
    allele_stats = compute_panel_allele_statistics(variant)
    panel_stats[var_id] = allele_stats

sys.stderr.write("Completed generating panel statistics.\n")
sys.stderr.write("Found %d variant IDs.\n"%(len(panel_stats)))
sys.stderr.write("\nReading Callset VCF.\n")
quals = [0,50,100,200]
for n, variant in enumerate(callset_reader):
    assert len(variant.ALT) == 1
    var_id = variant.INFO['ID']
    allele_stats, genotype_stats, counts = compute_callset_statistics(variant, qualities=quals, sample_pop_map=sample2pop, samples=callset_samples)
    assert list(allele_stats.keys()) == ['all', 'AFR', 'AMR', 'EAS', 'EUR', 'SAS']
    assert list(genotype_stats.keys()) == ['all', 'AFR', 'AMR', 'EAS', 'EUR', 'SAS']
    assert list(counts.keys()) == ['all', 'AFR', 'AMR', 'EAS', 'EUR', 'SAS']
    callset_stats[var_id] = [allele_stats, genotype_stats, counts]
    if (n+1)%1000 == 0:
        sys.stderr.write("\tRead %d variant records.\n"%(n))

sys.stderr.write("Completed generating callset statistics.\n")
sys.stderr.write("Qualities used for stat generation: %s.\n"%(','.join([str(q) for q in quals])))
sys.stderr.write("Found %d variant IDs.\n"%(len(callset_stats)))

# print stats for all IDs in genotypes VCF
header = [ 	'variant_id',
        'variant_type',
        'variant_length',
        'panel_allele_freq',
        'panel_alternative_alleles',
        'panel_total_alleles',
        'panel_unknown_alleles']

for c in ['all', 'AFR', 'AMR', 'EAS', 'EUR', 'SAS']:
    header.append('%s_allele_freq'%(c))
    header.append('%s_alternative_alleles'%(c))
    header.append('%s_total_alleles'%(c))
    header.append('%s_unknown_alleles'%(c))
    header.append('%s_heterozygosity'%(c))
    header.append('%s_heterozygous_genotypes'%(c))
    header.append('%s_homozygous_reference_genotypes'%(c))
    header.append('%s_homozygous_alternate_genotypes'%(c))
    header.append('%s_total_genotypes'%(c))

    for q in quals:
        header.append('%s_GQ>=%s'%(c,str(q)))

print('\t'.join(header))

assert len(panel_stats) == len(callset_stats)

for var_id in callset_stats:
    if not var_id in panel_stats:
        continue
    
    _,_,var_type,_,var_len = var_id.split('-')
    line = [var_id,
            var_type,
            var_len,
            panel_stats[var_id].af,
            panel_stats[var_id].ac,
            panel_stats[var_id].an,
            panel_stats[var_id].untyped]
    for c in ['all', 'AFR', 'AMR', 'EAS', 'EUR', 'SAS']:
        line.extend(
            [callset_stats[var_id][0][c].af,
            callset_stats[var_id][0][c].ac,
            callset_stats[var_id][0][c].an,
            callset_stats[var_id][0][c].untyped,
            callset_stats[var_id][1][c].heterozygosity,
            callset_stats[var_id][1][c].het,
            callset_stats[var_id][1][c].hom_ref,
            callset_stats[var_id][1][c].hom_alt,
            callset_stats[var_id][1][c].total])
        # add counts for GQs
        for q in quals:
            line.append(str(callset_stats[var_id][2][c][q]))
    assert len(line) == len(header)
    print('\t'.join(line))