import gzip
import argparse
import pysam

def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--debug", action="store_true", default=False, help="Print debug messages")
    parser.add_argument("--vcf", required=True, help="Input VCF.")
    parser.add_argument("--output", required=True, help="Output VCF.")
    parser.add_argument("--bgzip", action='store_true', help="Flag to bgzip the output VCF.")

    options = parser.parse_args()
    
    vcf = options.vcf
    reader = None
    if is_file_gzipped(vcf):
        reader = gzip.open(vcf, 'rt')
    else:
        reader = open(vcf, 'r')
    
    out = options.output
    writer = None
    if options.bgzip:
        writer = pysam.libcbgzf.BGZFile(out, 'wb')
    else:
        writer = open(out, "w")
    
    while True:
        line = reader.readline()
        if not line:
            break
        if line[0] == "#":
            write_to_file(line, writer)
            continue
        fields = line.rstrip().split('\t')
        chr = fields[0]
        pos = int(fields[1])
        if chr != "chrX":
            write_to_file(line, writer)
            continue
        if pos >= 10001 and pos <= 2781479:
            write_to_file(line, writer)
            continue
        if pos >= 155701383:
            write_to_file(line, writer)
            continue
    
def write_to_file(line, writer):
    try:
        writer.write(line)
    except TypeError:
        writer.write(str.encode(line))

def is_file_gzipped(src):
    with open(src, "rb") as inp:
        return inp.read(2) == b'\x1f\x8b'

if __name__ == "__main__":
    main()