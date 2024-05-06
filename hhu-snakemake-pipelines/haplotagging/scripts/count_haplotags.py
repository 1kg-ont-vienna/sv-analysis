'''
Counts the number of reads haplotaged and the total number of reads for the sample.
Reports summary of all samples and if --output is provided, lists the reads exclusively present in FASTA
'''

import argparse
import subprocess

def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", required=True, help="Comma-separated list of samples")
    parser.add_argument("--tsv", required=True, help="Path to the folder containing TSVs. The format is <Path to TSV>/<Sample Name>/<Sample Name>.tsv")
    parser.add_argument("--fasta", required=True, help="Path to the folder containing FASTAs.The format is <Path to FASTA>/<Sample Name>.fasta")
    parser.add_argument("--output", help="Output directory to store list of reads exclusively present in FASTA. (Summary in stdout)")
    
    options = parser.parse_args()

    samples = options.samples.split(',')

    print("#Sample\tReadsInFASTA\tReadsInTSV\tReadsOnlyInFASTA")
    for sample in samples:
        tsv_count, tsv_reads = extract_tsv_count(options.tsv, sample, options.output)
        fasta_count, fasta_reads = extract_fasta_count(options.fasta, sample, options.output)
        print("%s\t%d\t%d\t%d"%(sample, fasta_count, tsv_count, fasta_count-tsv_count))
        if options.output:
            intersection_reads = list(fasta_reads.difference(tsv_reads))
            out = options.output+"/%s.txt"%(sample)
            with open(out, 'w') as f:
                for read in intersection_reads:
                    f.write(read)
                    f.write("\n")


def extract_tsv_count(path, sample, output):
    if output == None:
        file = path+"/%s/%s.tsv"%(sample, sample)
        cmd = "wc -l %s"%(file)
        process = subprocess.run(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        num = int(process.stdout.rstrip().split()[0])-1
        return num, None
    file = path+"/%s/%s.tsv"%(sample, sample)
    count = 0
    reads = set()
    with open(file, 'r') as f:
        for line in f:
            if line[0] == "#":
                continue
            count += 1
            reads.add(line.rstrip().split()[0])
    return count, reads


def extract_fasta_count(path, sample, output):
    if output == None:
        file = path+"/%s.fasta"%(sample)
        cat_cmd = "cat %s"%(file)
        cat_process = subprocess.Popen(cat_cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        grep_cmd = "grep ^>"
        grep_process = subprocess.Popen(grep_cmd.split(), stdin=cat_process.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        wc_cmd = "wc -l"
        wc_process = subprocess.run(wc_cmd.split(), stdin=grep_process.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        num = int(wc_process.stdout.rstrip().split()[0])
        cat_process.stdout.close()
        grep_process.stdout.close()
        return num, None
    file = path+"/%s.fasta"%(sample)
    count = 0
    reads = set()
    with open(file, 'r') as f:
        for line in f:
            if line[0] != ">":
                continue
            count += 1
            reads.add(line.rstrip().split()[0][1:])
    return count, reads
    

if __name__ == '__main__':
    main()