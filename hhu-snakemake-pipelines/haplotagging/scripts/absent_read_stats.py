'''
Processes the list of read IDs which are not present in the TSV.
Finds out which category they belong to: unmapped, chromosomal, extra-chromosomal
Outputs the count summary and separate files for the three categories.
'''

import argparse
import pysam

def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--absent-reads", required=True, help="Path to list of absent reads")
    parser.add_argument("--bam", required=True, help="Path to BAM")
    parser.add_argument("--sample", required=True, help="Sample Name")
    parser.add_argument("--output", required=True, help="Output directory to store list of reads of different categories.")
    
    options = parser.parse_args()

    reads = extract_reads(options.absent_reads)

    chromosomes = ["chr%s"%(i) for i in [str(j) for j in range(1,23)]+['X']]
    reads_unmapped = []
    reads_extra = []
    bamfile = pysam.AlignmentFile(options.bam, 'rb')
    for alignment in bamfile.fetch(until_eof=True):
        # check for unmapped read
        if alignment.is_unmapped:
            assert alignment.query_name in reads
            reads_unmapped.append(alignment.query_name)
            reads.remove(alignment.query_name)
            continue
        # WhatsHap haplotagging ignores secondary and supplementary alignments for the haplotag writer
        if alignment.is_secondary or alignment.is_supplementary:
            continue
        # check for extra-chromosomal category
        if alignment.reference_name not in chromosomes:
            if alignment.query_name in reads:
                reads_extra.append(alignment.query_name)
                reads.remove(alignment.query_name)
    # Whatever reads are left in the list `reads` must have aligned to the chromosomes.

    # writing the summary stats
    # Category  Count
    writer = open("%s/%s_absent-read-stats.tsv"%(options.output,options.sample), 'w')
    writer.write("unmapped\t%d\n"%(len(reads_unmapped)))
    writer.write("extra-chromosomal\t%d\n"%(len(reads_extra)))
    writer.write("chromosomal\t%d\n"%(len(reads)))
    writer.close()

    # writing separate read lists
    writer = open("%s/%s_absent-read-unmapped.txt"%(options.output,options.sample), 'w')
    for read in reads_unmapped:
        writer.write(read)
        writer.write("\n")
    writer.close()
    writer = open("%s/%s_absent-read-extra-chromosomal.txt"%(options.output,options.sample), 'w')
    for read in reads_extra:
        writer.write(read)
        writer.write("\n")
    writer.close()
    writer = open("%s/%s_absent-read-chromosomal.txt"%(options.output,options.sample), 'w')
    for read in reads:
        writer.write(read)
        writer.write("\n")
    writer.close()



def extract_reads(path):
    reads=[]
    with open(path, 'r') as f:
        for read in f:
            reads.append(read.rstrip())
    return reads



if __name__ == '__main__':
    main()