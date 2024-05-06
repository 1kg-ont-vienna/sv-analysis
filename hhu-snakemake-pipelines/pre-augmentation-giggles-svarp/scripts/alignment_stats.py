import argparse
import logging
import resource
from collections import defaultdict
import gzip
import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)

def setup_logging(debug):
    handler = logging.StreamHandler()
    root = logging.getLogger()
    root.addHandler(handler)
    root.setLevel(logging.DEBUG if debug else logging.INFO)

def validate_arguments(options):
    pass

def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--debug", action="store_true", default=False, help="Print debug messages")
    parser.add_argument("--assembly-list", required=True, help="Text file with the list of assembly-to-graph GAF files. The GAF files should be sorted.")
    parser.add_argument("--output", required=True, help="Output directory to store all the stats and graphs.")
    
    options = parser.parse_args()
    setup_logging(options.debug)
    
    validate_arguments(options)

    gaf_files = []
    with open(options.assembly_list, 'r') as f:
        for line in f:
            gaf_files.append(line.rstrip())
    
    for gaf in gaf_files:
        get_stats(gaf, options.output)

    memory_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    logger.info("\nMemory Information")
    logger.info("  Maximum memory usage: %.3f GB", memory_kb / 1e6)
    

def get_stats(gaf, out):
    print("\nINFO: Working on file %s"%(gaf))
    if is_file_gzipped(gaf):
        reader = gzip.open(gaf, 'rt')
    else:
        reader = open(gaf, 'r')
    stats = defaultdict(lambda: [])
    total_alignments = 0
    total_alignment_len = 0
    num25 = 0
    num50 = 0
    num75 = 0
    num90 = 0
    num100 = 0
    num_prim = 0
    while True:
        alignment = reader.readline()
        if not alignment:
            break
        alignment = alignment.split('\t')
        seq_name = alignment[0]
        q_len = int(alignment[1])
        q_start = int(alignment[2])
        q_end = int(alignment[3])
        total_alignment_len += abs(q_start-q_end)
        percentage_aligned = abs(q_start-q_end)/q_len
        num_res_match = int(alignment[9])
        block_len = int(alignment[10])
        mapq = int(alignment[11])
        tags = defaultdict(lambda: None)
        for fields in alignment[12:]:
            if fields.startswith('tp:A:'):
                tags['tp'] = fields[5:]
            if fields.startswith('NM:i:'):
                tags['NM'] = int(fields[5:])
            if fields.startswith('dv:f:'):
                tags['dv'] = float(fields[5:])
        stats[seq_name].append([percentage_aligned, num_res_match, block_len, mapq, tags])
        total_alignments += 1
        if percentage_aligned < 0.25:
            num25 += 1
        elif percentage_aligned < 0.5:
            num50 += 1
        elif percentage_aligned < 0.75:
            num75 += 1
        elif percentage_aligned < 0.9:
            num90 += 1
        else:
            num100 += 1
        if tags['tp'] == 'P':
            num_prim += 1

    num_align = []
    # This stores the stat for the largest alignment per read
    per_read_stat = {'p_align': [], 'num_res_match': [], 'mapq': [], 'dv': []}
    total_stat = {'p_align': [], 'num_res_match': [], 'mapq': [], 'dv': []}
    for read, stat in stats.items():
        num_align.append(len(stat))
        max_len = 0
        index = None
        for n,s in enumerate(stat):
            total_stat['p_align'].append(s[0])
            total_stat['num_res_match'].append(s[1])
            total_stat['mapq'].append(s[3])
            total_stat['dv'].append(s[4]['dv'])
            if s[0] > max_len:
                index = n
        per_read_stat['p_align'].append(stat[index][0])
        per_read_stat['num_res_match'].append(stat[index][1])
        per_read_stat['mapq'].append(stat[index][3])
        per_read_stat['dv'].append(stat[index][4]['dv'])

    
    print("\tFound %d alignments (where %d are primary alignments) from %d sequences."%(total_alignments, num_prim, len(list(stats.keys()))))
    print("\tAlignment with percentage query length between:")
    print("\t\t0.00 and 0.25: %d"%(num25))
    print("\t\t0.25 and 0.50: %d"%(num50))
    print("\t\t0.50 and 0.75: %d"%(num75))
    print("\t\t0.75 and 0.90: %d"%(num90))
    print("\t\t0.90 and 1.00: %d"%(num100))
    print("\tTotal base pairs in alignment (not corrected): %d"%(total_alignment_len))

    graph(per_read_stat, total_stat, num_align, gaf, out)

def graph(read, total, align, gaf, out):
    
    hap = "_".join(gaf.split("/")[-1].split(".")[0:2])
    for y_label in read.keys():
        read_data = read[y_label]
        total_data = total[y_label]
        p = hap+"_"+y_label
        fig, axes = plt.subplots(2, 1, figsize=(15, 15))
        plt.suptitle(y_label)
        axes[0].hist(read_data, bins=20, rwidth = 0.8)
        axes[1].hist(total_data, bins=20, rwidth = 0.8)
        axes[0].set_title('Per Read Data (Selected Alignment with Highest Alignment Length)')
        axes[1].set_title('Data for all the alignments')
        plt.savefig(out+p+'.png')

    fig, axes = plt.subplots(1, 1, figsize=(15, 10))
    axes.hist(align, bins=20, rwidth = 0.8)
    axes.set_title('Number of Alignments per read')
    plt.savefig('%s%s_num_align.png'%(out,hap))


def is_file_gzipped(src):
    with open(src, "rb") as inp:
        return inp.read(2) == b'\x1f\x8b'

if __name__ == "__main__":
    main()