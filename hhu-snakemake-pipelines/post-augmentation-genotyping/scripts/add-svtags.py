import sys

for line in sys.stdin:
    if line.startswith('##'):
        print(line.strip())
        continue
    if line.startswith('#'):
        print("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of structural variation\">")
        print("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variation\">")
        print(line.strip())
        continue
    fields = line.strip().split()


    # determine SVLEN and SVTYPE
    ref_allele = fields[3]
    alt_alleles = fields[4].split(',')
    assert len(alt_alleles) == 1

    length = max([len(ref_allele), len(alt_alleles[0])])

    if len(ref_allele) < len(alt_alleles[0]):
        svtype = "INS"
    elif len(ref_allele) > len(alt_alleles[0]):
        svtype = "DEL"
    else:
        svtype = "OTHER"
    
    fields[7] = fields[7] + ';SVTYPE=' + svtype
    fields[7] = fields[7] + ';SVLEN=' + str(length)
    print('\t'.join(fields))