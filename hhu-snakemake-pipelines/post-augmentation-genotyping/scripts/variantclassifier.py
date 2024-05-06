from enum import Enum

class VariantType(Enum):
	snp = 0
	small_insertion = 1
	small_deletion = 2
	small_complex = 3
	midsize_insertion = 4
	midsize_deletion = 5
	midsize_complex = 6
	large_insertion_50_1000 = 7
	large_deletion_50_1000 = 8
	large_complex_50_1000 = 9
	large_insertion_1000_10000 = 10
	large_deletion_1000_10000 = 11
	large_complex_1000_10000 = 12
	large_insertion_10000_50000 = 13
	large_deletion_10000_50000 = 14
	large_complex_10000_50000 = 15
	large_insertion_50000 = 16
	large_deletion_50000 = 17
	large_complex_50000 = 18


def determine_class_from_alleles(ref_allele, alt_allele):
	if len(ref_allele) == len(alt_allele) == 1:
		return 'SNV'
	if (len(ref_allele) > 1) and (len(alt_allele) == 1):
		return 'DEL'
	if (len(ref_allele) == 1) and (len(alt_allele) > 1):
		return 'INS'
	return 'COMPLEX'


def determine_type_from_alleles(ref_allele, alt_allele):
	alleles = [ref_allele, alt_allele]

	if all([len(a) == 1 for a in alleles]):
		return VariantType.snp

	is_deletion = (len(alleles[0]) > 1) and (len(alleles[1]) == 1)
	is_insertion = (len(alleles[0]) == 1) and (len(alleles[1]) > 1)

	varlen = max([len(a) for a in alleles])
	if is_deletion or is_insertion:
		assert len(alleles) == 2
		varlen = varlen - 1

	if varlen < 20:
		if is_insertion:
			return VariantType.small_insertion
		if is_deletion:
			return VariantType.small_deletion
		return VariantType.small_complex

	if varlen >= 20 and varlen < 50:
		if is_insertion:
			return VariantType.midsize_insertion
		if is_deletion:
			return VariantType.midsize_deletion
		return VariantType.midsize_complex

	if varlen >= 50:
		if is_insertion:
			return VariantType.large_insertion
		if is_deletion:
			return VariantType.large_deletion
		return VariantType.large_complex


def determine_type_from_id(varid):
	if 'SNV' in varid:
		return VariantType.snp

	is_deletion = 'DEL' in varid
	is_insertion = 'INS' in varid

	varlen = int(varid.split('-')[-1])
	
	if varlen < 20:
		if is_insertion:
			return VariantType.small_insertion
		if is_deletion:
			return VariantType.small_deletion
		return VariantType.small_complex

	if varlen >= 20 and varlen < 50:
		if is_insertion:
			return VariantType.midsize_insertion
		if is_deletion:
			return VariantType.midsize_deletion
		return VariantType.midsize_complex

	if varlen >= 50 and varlen < 1000:
		if is_insertion:
			return VariantType.large_insertion_50_1000
		if is_deletion:
			return VariantType.large_deletion_50_1000
		return VariantType.large_complex_50_1000
	
	if varlen >= 1000 and varlen < 10000:
		if is_insertion:
			return VariantType.large_insertion_1000_10000
		if is_deletion:
			return VariantType.large_deletion_1000_10000
		return VariantType.large_complex_1000_10000

	if varlen >= 10000 and varlen < 50000:
		if is_insertion:
			return VariantType.large_insertion_10000_50000
		if is_deletion:
			return VariantType.large_deletion_10000_50000
		return VariantType.large_complex_10000_50000

	if varlen >= 50000:
		if is_insertion:
			return VariantType.large_insertion_50000
		if is_deletion:
			return VariantType.large_deletion_50000
		return VariantType.large_complex_50000


def determine_variant_type(line):
	"""
	Determines the variant type based on
	the IDs of the variant alleles.
	"""
	fields = line.split()
	info_fields = {f.split('=')[0] : f.split('=')[1] for f in fields[7].split(';') if '=' in f}
	assert 'ID' in info_fields
	allele_ids = info_fields['ID'].split(',')
	# handle merged IDs
	all_ids = []
	for i in allele_ids:
		for j in i.split(':'):
			all_ids.append(j)
	all_ids = list(set(all_ids))
	assert len(all_ids) > 0
	
	if all(['SNV' in i for i in all_ids]) and all([ i.split('-')[1] == all_ids[0].split('-')[1]  for i in all_ids]):
		# all SNVs starting at the same base
		return VariantType.snp

	is_deletion = (len(all_ids) == 1) and 'DEL' in all_ids[0]
	is_insertion = (len(all_ids) == 1) and 'INS' in all_ids[0]

	alleles = [fields[3]] + [f for f in fields[4].split(',')]
	varlen = max([len(a) for a in alleles])
	if is_deletion or is_insertion:
		assert len(all_ids) == 1
		varlen = int(all_ids[0].split('-')[-1])

	if varlen < 20:
		if is_insertion:
			return VariantType.small_insertion
		if is_deletion:
			return VariantType.small_deletion
		return VariantType.small_complex

	if varlen >= 20 and varlen < 50:
		if is_insertion:
			return VariantType.midsize_insertion
		if is_deletion:
			return VariantType.midsize_deletion
		return VariantType.midsize_complex

	if varlen >= 50 and varlen < 1000:
		if is_insertion:
			return VariantType.large_insertion_50_1000
		if is_deletion:
			return VariantType.large_deletion_50_1000
		return VariantType.large_complex_50_1000
	
	if varlen >= 1000 and varlen < 10000:
		if is_insertion:
			return VariantType.large_insertion_1000_10000
		if is_deletion:
			return VariantType.large_deletion_1000_10000
		return VariantType.large_complex_1000_10000

	if varlen >= 10000 and varlen < 50000:
		if is_insertion:
			return VariantType.large_insertion_10000_50000
		if is_deletion:
			return VariantType.large_deletion_10000_50000
		return VariantType.large_complex_10000_50000

	if varlen >= 50000:
		if is_insertion:
			return VariantType.large_insertion_50000
		if is_deletion:
			return VariantType.large_deletion_50000
		return VariantType.large_complex_50000


def determine_variant_from_line(line):
	"""
	Determines the variant type directly from the alleles.
	"""
	splitted = line.split()
	alleles = [splitted[3]] + [s for s in splitted[4].split(',')]

	if all([len(a) == 1 for a in alleles]):
		return VariantType.snp

	is_deletion = (len(alleles[0]) > 1) and (len(alleles[1]) == 1) and (len(alleles) < 3)
	is_insertion = (len(alleles[0]) == 1) and (len(alleles[1]) > 1) and (len(alleles) < 3)

	varlen = max([len(a) for a in alleles])
	if is_deletion or is_insertion:
		assert len(alleles) == 2
		varlen = varlen - 1

	if varlen < 20:
		if is_insertion:
			return VariantType.small_insertion
		if is_deletion:
			return VariantType.small_deletion
		return VariantType.small_complex

	if varlen >= 20 and varlen < 50:
		if is_insertion:
			return VariantType.midsize_insertion
		if is_deletion:
			return VariantType.midsize_deletion
		return VariantType.midsize_complex

	if varlen >= 50 and varlen < 1000:
		if is_insertion:
			return VariantType.large_insertion_50_1000
		if is_deletion:
			return VariantType.large_deletion_50_1000
		return VariantType.large_complex_50_1000
	
	if varlen >= 1000 and varlen < 10000:
		if is_insertion:
			return VariantType.large_insertion_1000_10000
		if is_deletion:
			return VariantType.large_deletion_1000_10000
		return VariantType.large_complex_1000_10000

	if varlen >= 10000 and varlen < 50000:
		if is_insertion:
			return VariantType.large_insertion_10000_50000
		if is_deletion:
			return VariantType.large_deletion_10000_50000
		return VariantType.large_complex_10000_50000

	if varlen >= 50000:
		if is_insertion:
			return VariantType.large_insertion_50000
		if is_deletion:
			return VariantType.large_deletion_50000
		return VariantType.large_complex_50000


def determine_pyvcf_type(record):
	"""
	Determine variant type from pyvcf
	record which does not have IDs.
	"""
	alleles = [record.REF] + record.ALT

	if record.is_snp:
		return VariantType.snp

	is_deletion = (len(record.ALT) == 1) and (len(record.REF) > 1) and (len(record.ALT[0]) == 1)
	is_insertion = (len(record.ALT) == 1) and (len(record.REF) == 1) and (len(record.ALT[0]) > 1)
	
	varlen = max([len(a) for a in alleles])
	if is_deletion or is_insertion:
		# insertions and deletions always have one base more
		varlen -= 1

	assert varlen >= 1

	if varlen < 20:
		if is_insertion:
			return VariantType.small_insertion
		if is_deletion:
			return VariantType.small_deletion
		return VariantType.small_complex

	if varlen >= 20 and varlen < 50:
		if is_insertion:
			return VariantType.midsize_insertion
		if is_deletion:
			return VariantType.midsize_deletion
		return VariantType.midsize_complex

	if varlen >= 50 and varlen < 1000:
		if is_insertion:
			return VariantType.large_insertion_50_1000
		if is_deletion:
			return VariantType.large_deletion_50_1000
		return VariantType.large_complex_50_1000
	
	if varlen >= 1000 and varlen < 10000:
		if is_insertion:
			return VariantType.large_insertion_1000_10000
		if is_deletion:
			return VariantType.large_deletion_1000_10000
		return VariantType.large_complex_1000_10000

	if varlen >= 10000 and varlen < 50000:
		if is_insertion:
			return VariantType.large_insertion_10000_50000
		if is_deletion:
			return VariantType.large_deletion_10000_50000
		return VariantType.large_complex_10000_50000

	if varlen >= 50000:
		if is_insertion:
			return VariantType.large_insertion_50000
		if is_deletion:
			return VariantType.large_deletion_50000
		return VariantType.large_complex_50000


