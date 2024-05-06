import sys, argparse
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

def plot_distances(all_variants, vartypes, outname):
    """
        Plots the distances between variants of the same
        type.
    """
    with PdfPages(outname) as pdf:
        for vartype in vartypes:
            variants = [v for v in all_variants if v.vartype == vartype]
            sys.stderr.write('Plotting histogram for ' + str(len(variants)) + ' ' + vartype + '...\n')
            plt.figure()
            plt.title('variant type: ' + vartype)
            distances = []
            for i in range(0, len(variants)-1):
                if variants[i].chrom != variants[i+1].chrom:
                    continue
                pos_a = variants[i].start
                pos_b = variants[i+1].start
                assert pos_a <= pos_b
                distances.append(pos_b - pos_a)
            plt.hist(distances, bins=50, range=(0,500))
            plt.ylabel('Count')
            plt.yscale('log')
            plt.xlabel('Distance to next variant (bp)');
            pdf.savefig()
            plt.close()
        

def record_from_cluster(cluster, callsets):
    """
        Given variants that all represent the same site,
        construct csv record. Coordinates are defined by
        callset with highest confidence indicated by
        the order of elements in callsets.
    """
    var_id = None
    chrom = None
    start = None
    end = None
    length = None
    vartype = None
    qual = None
    callset_to_variant = {c:[] for c in callsets}
#    assert all([c.vartype == cluster[0].vartype for c in cluster])
    for c in cluster:
        callset_to_variant[c.callset].append(c)

    for c in callsets:
        if callset_to_variant[c] != []:
            var_id = callset_to_variant[c][0].id
            chrom = callset_to_variant[c][0].chrom
            start = callset_to_variant[c][0].start
            end = callset_to_variant[c][0].end
            length = callset_to_variant[c][0].varlen
            vartype = callset_to_variant[c][0].vartype
            qual = callset_to_variant[c][0].quality
            read_depth = callset_to_variant[c][0].read_depth
            allele_depth = callset_to_variant[c][0].allele_depth
            break

    callset_to_id = ['nan'] * len(callsets)
    callset_to_present = ['False'] * len(callsets)
    for i, c in enumerate(callsets):
        if callset_to_variant[c] != []:
            callset_to_id[i] = ';'.join([v.id for v in callset_to_variant[c]])
            callset_to_present[i] = 'True'
    tab_record = [var_id, chrom, str(start), str(end), str(length), vartype, qual, read_depth, allele_depth] + callset_to_present + callset_to_id
    vcf_record = [chrom, str(start), var_id, 'N', '<' + vartype + '>', qual, 'PASS', 'END=' + str(end) + ';SVTYPE=' + str(vartype) + ';SVLEN=' + str(length)]

    return tab_record, vcf_record


class Variant:

    def __init__(self, chrom, start, end, varlen, vartype, callset, quality, read_depth='nan', allele_depth='nan'):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.varlen = varlen
        self.vartype = vartype
        self.callset = callset
        self.quality = str(quality)
        self.read_depth = str(read_depth)
        self.allele_depth = str(allele_depth)
        self.id = '-'.join([self.chrom, str(self.start), self.vartype, str(self.varlen), self.callset])

    def __lt__(self, other):
        if self.chrom != other.chrom:
           return self.chrom < other.chrom
        else:
           return self.start < other.start

    def __len__(self):
        return self.end - self.start

    def __str__(self):
        return ','.join(['Variant(' + str(self.chrom), str(self.start), str(self.end), str(self.varlen), self.vartype, self.callset, self.quality + ')'])

    def __eq__(self, other):
        if self.chrom != other.chrom:
            return False
        if self.start != other.start:
            return False
        if self.end != other.end:
            return False
        if self.varlen != other.varlen:
            return False
        if self.vartype != other.vartype:
            return False
        if self.callset != other.callset:
            return False
        if self.quality != other.quality:
            return False
        if self.read_depth != other.read_depth:
            return False
        if self.allele_depth != other.allele_depth:
            return False
        return True

    def __hash__(self):
        return hash((self.chrom, self.start, self.end, self.varlen, self.vartype, self.callset, self.quality. self.read_depth, self.allele_depth))


def reciprocal_overlap(var1, var2, r):
    """
        Checks whether both variants have
        a reciprocal overlap of at least r.
    """
    if var1.start > var2.end or var2.start > var1.end:
        return False
    overlap = min(var1.end, var2.end) - max(var1.start, var2.start)
    if overlap > r * len(var1):
        if overlap > r * len(var2):
            return True
    return False


def coordinates_match(var1, var2):
    """
        Coordinates match if deviations of start and
        end are at most 200 bp and the length differs
        by at most 25%.
    """
    start_distance = abs(var1.start - var2.start)
    end_distance = abs(var1.end - var2.end)
    size_diff = var2.varlen / var1.varlen
    
    if (start_distance < 200) and (end_distance < 200) and (0.5 < size_diff < 1.5):
        return True
    else:
        return False


def same_variant(var1, var2, only_different_callsets=False):
    """
        Decides whether two variants are the same
        and can be merged.
    """
    if only_different_callsets and (var1.callset == var2.callset):
        return False
    assert var1.chrom == var2.chrom
    # variants need to be of same type
    if var1.vartype != var2.vartype:
        return False
    if reciprocal_overlap(var1, var2, 0.5) or coordinates_match(var1, var2):
        return True
    else:
        return False


def merge_variants(variants, callset_names, only_different_callsets = False):
    """
        Given a set of overlapping variants, merge them
        in the order of the given callset_names
    """
    callset_to_variants = defaultdict(list)
    for v in variants:
        callset_to_variants[v.callset].append(v)

    ordered_variants = []
    for callset in callset_names:
        if callset_to_variants[callset] != []:
            ordered_variants.extend(callset_to_variants[callset])

    clusters = [[ordered_variants[0]]]
    for v in ordered_variants[1:]:
        overlap_exists = False
        for c in range(len(clusters)):
            if same_variant(v, clusters[c][0], only_different_callsets):
                overlap_exists = True
                clusters[c].append(v)

        if not overlap_exists:
            clusters.append([v])

    tsv_results = []
    vcf_results = []
    for cluster in clusters:
        tsv_result, vcf_result = record_from_cluster(cluster, callset_names)
        tsv_results.append(tsv_result)
        vcf_results.append(vcf_result)
    return tsv_results, vcf_results



def find_clusters(variants, offset=0):
    """
        Finds sets of overlapping variants and
        merges such that are the same.
    """
    assert len(variants) > 1
    current_start = variants[0].start
    current_end = variants[0].end
    current_cluster = [variants[0]]
    current_chrom = variants[0].chrom
    for v in variants[1:]:
        if (v.start > (current_end + offset)) or (v.chrom != current_chrom):
            yield current_cluster
            current_cluster = []
        current_chrom = v.chrom
        current_cluster.append(v)
        current_start = v.start
        current_end = max(v.end, current_end)

    if len(current_cluster) != 0:
        yield current_cluster


def parse_vcf(filename, name, variants, vartypes, id_from_vcf=False):
    """
        Reads a VCF and stores Variant objects
        in variants.
    """
    var_counter = 0
    for line in open(filename, 'r'):
        if line.startswith('#'):
            continue
        fields = line.split()
        chrom = fields[0]
        start = int(fields[1])
        quality = fields[5]
        info_field = {i.split('=')[0]:i.split('=')[1] for i in fields[7].split(';') if '=' in i if '=' in i}
        assert 'SVTYPE' in info_field
        vartype = info_field['SVTYPE']
        if vartype == "BND":
            # translocations cannot be handled and will be skipped
            sys.stderr.write('Skipping variant at position ' + chrom + ':' + str(start) + ' in ' + filename + ' since variant type is BND.\n')
            continue
        if not 'END' in info_field:
            # make sure the VCF is not symbolic in this case
            assert (not '<' in fields[4])
            end = start + len(fields[3])
            assert end >= start
        else:
            end = int(info_field['END'])
        if not 'SVLEN' in info_field:
            if vartype == 'INV':
                varlen = abs(end - start)
            elif vartype == 'INS':
                assert 'SVINSLEN' in info_field
                varlen = abs(int(info_field['SVINSLEN']))
            else:            
                sys.stderr.write('Skipping variant at position ' + chrom + ':' + str(start) + ' in ' + filename + ' since no SVLEN given.\n')
                continue
        else:
            varlen = abs(int(info_field['SVLEN']))

        if varlen < 50:
            sys.stderr.write('Skipping variant at position ' + chrom + ':' + str(start) + ' in ' + filename + ' since its length is smaller than 50bp (no an SV).\n')
            continue
        gts = {k : v for k,v in zip(fields[8].split(':'), fields[9].split(':'))}

        # determine read depth. Depending on caller, this is encoded in different fields
        read_depth = 'nan'
        allele_depth = 'nan'
        if ('DP' in gts) or ('DR' in gts and 'DV' in gts):
            if 'DP' in gts:
                 read_depth = gts['DP']
            else:
                assert 'DR' in gts
                assert 'DV' in gts
                if gts['DR'] != '.' and gts['DV'] != '.':
                    read_depth = str(int(gts['DR']) + int(gts['DV']))

            # determine allele depth

            if 'DV' in gts:
                allele_depth = gts['DV']
            else:
                assert 'AD' in gts
                ad = gts['AD'].split(',')
                assert len(ad) == 2
                allele_depth = ad[1]

        variant = Variant(chrom, start, end, varlen, vartype, name, quality, read_depth, allele_depth)
        if id_from_vcf:
            variant.id = fields[2]
        variants.append(variant)
        vartypes.add(vartype)
        var_counter += 1

    return var_counter


def parse_table(filename, variants, fulltable):
    """
        Reads a table constructed by intersect
        and stores variant objects.
    """
    var_counter = 0
    for line in open(filename, 'r'):
        if line.startswith('ID'):
            fulltable['header'] = line.strip()
            continue
        fields = line.split()
        chrom = fields[1]
        start = int(fields[2])
        end = int(fields[3])
        varlen = int(fields[4])
        vartype = fields[5]
        quality = fields[6]
        read_depth = fields[7]
        allele_depth = fields[8]
        variant = Variant(chrom, start, end, varlen, vartype, 'table', quality, read_depth, allele_depth)
        # set ID to the one given in table
        variant.id = fields[0] + '-table'
        variants.append(variant)
        fulltable[fields[0] + '-table'] = fields
        var_counter += 1
    return var_counter


def run_intersect(callsets, names, outtable, outvcf, outpdf):
    """
        Runs intersect command.
    """   
    sys.stderr.write('callset merging order: ' + '>'.join(names) + '\n')
    variants = []
    vartypes = set([])
    for callset, name in zip(callsets, names):
        count = parse_vcf(callset, name, variants, vartypes)
        sys.stderr.write('Read ' + str(count) + ' variants from ' + callset + '.\n')
    variants.sort()
    # plot distances between variants of same type
    plot_distances(variants, vartypes, outpdf)
    with open(outtable, 'w') as outtable, open(outvcf, 'w') as outvcf:
        tsv_header = [ 'ID', 'chromosome', 'start', 'end', 'length', 'var_type', 'quality', 'read_depth', 'allele_depth'] + ['in_' + c for c in names] + ['ID_' + c for c in names]
        outtable.write('\t'.join(tsv_header) + '\n')
        outvcf.write("##fileformat=VCFv4.2\n")
        outvcf.write("##ALT=<ID=DEL,Description=\"Deletion\">\n")
        outvcf.write("##ALT=<ID=DUP,Description=\"Duplication\">\n")
        outvcf.write("##ALT=<ID=INV,Description=\"Inversion\">)\n")
        outvcf.write("##ALT=<ID=INVDUP,Description=\"InvertedDUP with unknown boundaries\">\n")
        outvcf.write("##ALT=<ID=TRA,Description=\"Translocation\">ņ")
        outvcf.write("##ALT=<ID=INS,Description=\"Insertion\">\n")
        outvcf.write("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of the SV\">\n")
        outvcf.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n")
        outvcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        nr_records = 0
        for cluster in find_clusters(variants, 200):
            tab_records, vcf_records = merge_variants(cluster, names, True)
            for r,v in zip(tab_records, vcf_records):
                outtable.write('\t'.join(r) + '\n')
                outvcf.write('\t'.join(v) + '\n')
                nr_records += 1
    sys.stderr.write('Number of variants in union: ' + str(nr_records) + '.\n')


def run_annotate(table, callset, name):
    sys.stderr.write('Adding callset ' + callset + ' to table ' + table + '.\n')
    variants = []
    vartypes = set([])
    fulltable = {}
    count_table = parse_table(table, variants, fulltable)
    sys.stderr.write('Read ' + str(count_table) + ' variants from ' + table + '.\n')
    
    count_callset = parse_vcf(callset, name, variants, vartypes, True)
    sys.stderr.write('Read ' + str(count_callset) + ' variants from ' + callset + '.\n')

    variants.sort()
    names = ['table', name]
    nr_records = 0
    # print header
    print(fulltable['header'] + '\tin_' + name + '\tID_' + name)
    for cluster in find_clusters(variants, 200):
        tab_records, vcf_records = merge_variants(cluster, names, True)
        for r in tab_records:
            if not 'table' in r[0]:
                continue
            assert r[0] in fulltable
            record = fulltable[r[0]] + [r[-3], r[-1]]
            print('\t'.join(record))
            nr_records += 1
    sys.stderr.write('Number of records: ' + str(nr_records) + '.\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='intersect-callsets.py', description=__doc__)
    subparsers = parser.add_subparsers(dest="subparser_name")

    parser_intersect = subparsers.add_parser('intersect', help='Intersect a list of VCF callsets.')
    parser_intersect.add_argument('-c', '--callsets', nargs='+', required=True, default=[], help='Callsets in VCF format. VCFs must be symbolic. Order corresponds to merging order.')
    parser_intersect.add_argument('-t', '--table', required=True, help='Name of the output table.')
    parser_intersect.add_argument('-v', '--vcf', required=True, help='Name of the output files.')
    parser_intersect.add_argument('-p', '--pdf', required=True, help='Name of the output files.')
    parser_intersect.add_argument('-n', '--names', nargs='+', required=True, default=[], help='Callset names. One per given VCF file.')

    parser_annotate = subparsers.add_parser('annotate', help='Add column indicating whether variant was seen in a given callset.')
    parser_annotate.add_argument('-t', '--table', required=True, help='Table created by intersect.')
    parser_annotate.add_argument('-v', '--vcf', required=True, help='Callset in VCF format.')
    parser_annotate.add_argument('-n', '--name', required=True, help='Callset name.')
    args = parser.parse_args()
    
    if args.subparser_name == "intersect":
        run_intersect(args.callsets, args.names, args.table, args.vcf, args.pdf)
    elif args.subparser_name == "annotate":
        run_annotate(args.table, args.vcf, args.name)
