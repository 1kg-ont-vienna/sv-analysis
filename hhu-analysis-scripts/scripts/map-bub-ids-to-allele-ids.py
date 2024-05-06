import argparse
from cyvcf2 import VCF

parser = argparse.ArgumentParser(prog='map-bub-ids-to-allele-ids.py', description="Create a map between bubble IDs and allele IDs.")
parser.add_argument('-vcf', metavar='vcf', help='Multiallelic VCF (with allele decomposition)')
args = parser.parse_args()

vcf_reader = VCF(args.vcf)
print("#Bubble-ID\tAllele-IDs")
for record in vcf_reader:
    bub_id = str(record.ID)
    id_field = record.INFO['ID']
    allele_ids = set()
    for i in id_field.split(','):
        for j in i.split(':'):
            allele_ids.add(j)
    print("%s\t%s"%(bub_id, ','.join(list(allele_ids))))
    
vcf_reader.close()