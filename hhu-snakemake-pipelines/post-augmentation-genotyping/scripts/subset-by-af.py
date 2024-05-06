import argparse
import gzip

def get_af(line):
    fields = line.strip().split('\t')[9:]
    gts = [x.split(':')[0] for x in fields]
    ac = 0
    tot = 0
    for gt in gts:
        if '/' in gt:
            gt = gt.split('/')
            assert len(gt) == 2
            if gt[0] == '1':
                ac += 1
                tot += 1
            elif gt[0] == '0':
                tot += 1
            else:
                assert gt[0] == '.'
            if gt[1] == '1':
                ac += 1
                tot += 1
            elif gt[1] == '0':
                tot += 1
            else:
                assert gt[1] == '.'
        elif '|' in gt:
            gt = gt.split('|')
            assert len(gt) == 2
            if gt[0] == '1':
                ac += 1
                tot += 1
            elif gt[0] == '0':
                tot += 1
            else:
                assert gt[0] == '.'
            if gt[1] == '1':
                ac += 1
                tot += 1
            elif gt[1] == '0':
                tot += 1
            else:
                assert gt[1] == '.'
        else:
            assert gt in ['.', '0', '1']
            if gt == '1':
                ac += 1
                tot += 1
            elif gt == '0':
                tot += 1
    
    return ac/tot

parser = argparse.ArgumentParser(prog='subset-by-af.py', description="subset vcfs based on allele frequency of one of the vcf")
parser.add_argument('-panel', metavar='panel', help='Biallelic panel VCF (gzipped).')
parser.add_argument('-callset', metavar='callset', help='Biallelic callset VCF (gzipped).')
parser.add_argument('-mode', metavar='mode', help='"panel" for panel af, "callset" for callset af, and "both" for both af.')
parser.add_argument('-min-af', metavar='min_af', help='Minimum af.')
parser.add_argument('-max-af', metavar='max_af', help='Maximum af.')
parser.add_argument('-output-panel', metavar='output_panel', help='Output panel VCF.')
parser.add_argument('-output-callset', metavar='output_callset', help='Output callset VCF.')
args = parser.parse_args()

panel = args.panel
callset = args.callset
mode = args.mode
if mode not in ['panel', 'callset', 'both']:
    raise RuntimeError('Wrong Mode specified.')
min_af = float(args.min_af)
max_af = float(args.max_af)
if min_af < 0 or min_af > 1 or max_af < 0 or max_af > 1 or min_af > max_af:
    raise RuntimeError('Invalid allele frequency values.')
output_panel = args.output_panel
output_callset = args.output_callset

panel_writer = open(output_panel, 'w')
callset_writer = open(output_callset, 'w')

# Reading the headers and writing into the output file
panel_reader = gzip.open(panel, 'rt')
for line in panel_reader:
    line = line.strip()
    if line[0] == "#":
        print(line, file=panel_writer)
    else:
        break
panel_reader.close()
callset_reader = gzip.open(callset, 'rt')
for line in callset_reader:
    line = line.strip()
    if line[0] == "#":
        print(line, file=callset_writer)
    else:
        break
callset_reader.close()

# Reinitializing files
panel_reader = gzip.open(panel, 'rt')
callset_reader = gzip.open(callset, 'rt')

# Skipping lines until the records start
for line in panel_reader:
    if line[0:2] == "#C":
        break
for line in callset_reader:
    if line[0:2] == "#C":
        break

while True:
    panel_line = panel_reader.readline()
    callset_line = callset_reader.readline()
    if not panel_line:
        assert not callset_line
        break
    assert panel_line[0] != "#"
    assert callset_line[0] != "#"
    panel_fields = panel_line.split('\t')
    callset_fields = callset_line.split('\t')
    # assert that both have same record info
    assert panel_fields[0] == callset_fields[0]
    assert panel_fields[1] == callset_fields[1]
    is_pass = False
    if mode == 'panel':
        af = get_af(panel_line.strip())
        if af >= min_af and af <= max_af:
            is_pass = True
    elif mode == 'callset':
        af = get_af(callset_line.strip())
        if af >= min_af and af <= max_af:
            is_pass = True
    elif mode == 'both':
        af1 = get_af(panel_line.strip())
        af2 = get_af(callset_line.strip())
        if af1 >= min_af and af1 <= max_af and af2 >= min_af and af2 <= max_af:
            is_pass = True
    if is_pass:
        print(panel_line.strip(), file=panel_writer)
        print(callset_line.strip(), file=callset_writer)

panel_writer.close()
callset_writer.close()