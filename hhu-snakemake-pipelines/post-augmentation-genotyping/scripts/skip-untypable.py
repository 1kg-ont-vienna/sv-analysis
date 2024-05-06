import sys
from collections import defaultdict

# list of IDs to skip
var_ids_file = sys.argv[1]
var_ids = defaultdict(lambda: False)

for line in open(var_ids_file, 'r'):
	var_ids[line.strip()] = True

sys.stderr.write('Read ' + str(len(var_ids)) + ' IDs from input list.\n')

total_lines = 0
skipped_lines = 0
for line in sys.stdin:
	if line.startswith('#'):
		print(line[:-1])
		continue
	fields = line.split()
	info_fields = { i.split('=')[0] : i.split('=')[1].strip() for i in fields[7].split(';') if '=' in i}
	assert 'ID' in info_fields
	assert len(info_fields['ID'].split(',')) == 1
	var_id = info_fields['ID']
	total_lines += 1
	if var_ids[var_id]:
		skipped_lines += 1
		continue
	print(line[:-1])

sys.stderr.write('Skipped ' + str(skipped_lines) + ' of ' + str(total_lines) + ' from input vcf.\n')
