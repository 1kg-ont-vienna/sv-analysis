# cat *.log | grep -A11 ": all" | cut -d '     ' -f 2-5 | python ~/Work/1000GP/genotyping-results-localcopy/scripts/mendel-stats-arrange.py

import sys

out = {x: '' for x in range(4,13)}
for n, line in enumerate(sys.stdin):
    line = line.rstrip()
    if n%13 > 3 and n%13 < 12:
        out[n%13] += line+'\t'

for x in out.keys():
    print(out[x])