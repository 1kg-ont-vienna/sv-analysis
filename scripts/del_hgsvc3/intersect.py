#! /usr/bin/env python

from __future__ import print_function
import csv
import argparse
import collections

parser = argparse.ArgumentParser(description='SV intersection')
parser.add_argument('-a', '--fileA', metavar='a.bed', required=True, dest='fileA', help='file A')
parser.add_argument('-b', '--fileB', metavar='b.bed', required=True, dest='fileB', help='file B')
parser.add_argument('-v', '--unmatched', dest='unmatched', action='store_true', default=False, help='unmatched SVs')
args = parser.parse_args()

win = 500
sizeRatio = 0.25
recOverlap = 0.25
df = collections.defaultdict(list)
with open(args.fileB) as f:
    reader = csv.reader(f, delimiter="\t")
    for row in reader:
        df[row[0]].append({ 'start': int(row[1]), 'end': int(row[2]), 'id': row[3] })


with open(args.fileA) as f:
    reader = csv.reader(f, delimiter="\t")
    for row in reader:
        start = int(row[1])
        end = int(row[2])
        l1 = end - start
        matched = False
        for comp in df[row[0]]:
            if (comp['end'] + win < start) or (end + win < comp['start']):
                continue
            l2 = comp['end'] - comp['start']
            if l1 > l2:
                tmp = l1
                l1 = l2
                l2 = tmp
            if (abs(start - comp['start']) >= win) or (abs(end - comp['end']) >= win) or (l1 / l2 < sizeRatio):
                # Check reciprocal overlap
                if (comp['end'] < start) or (end < comp['start']):
                    continue
                (ovS, ovE) = sorted([start, end, comp['start'], comp['end']])[1:3]
                iSize = ovE - ovS
                if (iSize/l1 <= recOverlap) or (iSize/l2 <= recOverlap):
                    continue
            matched = True
            if not args.unmatched:
                print('\t'.join(row), row[0], str(comp['start']), str(comp['end']) , comp['id'], sep="\t")
        if args.unmatched:
            if not matched:
                print('\t'.join(row), sep="\t")
