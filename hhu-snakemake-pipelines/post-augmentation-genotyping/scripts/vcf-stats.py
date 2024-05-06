import sys

nr_snps = 0
nr_indels = 0
nr_svs = 0

for line in sys.stdin:
    if line.startswith('#'):
        continue
    fields = line.strip().split()
    ref_allele = fields[3]
    alt_alleles = fields[4].split(',')
    if len(alt_alleles) > 1:
        raise RuntimeError("VCF file must be biallelic.")
    assert len(alt_alleles) == 1
    # determine variant class
    if (len(ref_allele) == 1) and (len(alt_alleles[0]) == 1):
        nr_snps += 1
        continue

    if (len(ref_allele) < 50) and (len(alt_alleles[0]) < 50):
        nr_indels += 1
        continue

    nr_svs += 1

print('number of SNPs: ' + str(nr_snps))
print('number of indels (< 50 bp): ' + str(nr_indels))
print('number of SVs (>= 50 bp): ' + str(nr_svs))   