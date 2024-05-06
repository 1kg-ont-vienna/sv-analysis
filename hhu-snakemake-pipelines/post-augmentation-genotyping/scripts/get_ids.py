import sys

for line in sys.stdin:
	if line.startswith('#'):
		continue
	fields = line.split()
	info_fields = {i.split('=')[0] : i.split('=')[1] for i in fields[7].split(';') if '=' in i}
	assert 'ID' in info_fields
	assert len(fields[4].split(',')) == 1
	print('-'.join([fields[0], fields[1], fields[3], fields[4]]))
