import sys

for line in sys.stdin:
    if line.startswith('#'):
        print(line.strip())
        continue
    fields = line.split()
    if fields[0] == 'chrX' or fields[0] == 'chrY':
        continue
    if fields[6] == '.':
        fields[6] = 'PASS'
    print('\t'.join(fields))