import json
import csv
import argparse

parser = argparse.ArgumentParser(description='Insertion reads')
parser.add_argument('-w', '--window', metavar='100', required=False, dest='window', help='window size')
parser.add_argument('-c', '--csv', metavar='test.csv', required=True, dest='csv', help='input csv file')
args = parser.parse_args()

win = 100
if args.window:
    win = int(args.window)

print("chrom", "start", "end", "type", "sample", sep="\t")
with open(args.csv) as f:
    reader = csv.DictReader(f)
    for row in reader:
        for k in row.keys():
            if (k.startswith("HG")) or (k.startswith("NA")):
                gt = 0
                for val in row[k].split('|'):
                    if val == '1':
                        gt += 1
                if gt > 0:
                    startval = int(row['POS']) - win
                    endval = int(row['POS']) + win
                    if (row["CHROM"] != "chrX") and (row["CHROM"] != "chrY"):
                        print(row["CHROM"], startval, endval, row['TE_Designation'], k, sep='\t')
