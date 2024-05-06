import pysam
import math

new_vcf = ""

class VCF:
    def __init__(self, chr_name, start_pos, var_id, ref=None, alt=None, qual=None, filt=None, svtype=None, svlen=None, contig_id=None, genotype=None):
        self.chr_name = chr_name
        self.start_pos = start_pos
        self.end_pos = start_pos + svlen
        self.var_id = var_id
        self.ref = ref
        self.alt = alt
        self.qual = qual
        self.filt = filt
        self.svtype = svtype
        self.svlen = svlen
        self.contig_id = contig_id
        self.genotype = genotype
        self.hasDups = False


class Contig:
    def __init__(self, chr_name, ref=None, alt=None, genotype=None):
        self.chr_name = chr_name
        self.start_pos = math.inf
        self.end_pos = -math.inf
        self.ref = ref
        self.alt = alt
        self.genotype = genotype
        self.vars = []


def lists_overlap(a, b):
  sb = set(b)
  return any(el in sb for el in a)


def is_file_gzipped(src):
    with open(src, "rb") as inp:
        return inp.read(2) == b'\x1f\x8b'


def is_overlapping(x_start, x_end, y_start, y_end):
    return max(x_start, y_start) <= min(x_end, y_end)


def parse_vcf(filename):
    import re
    import gzip  
    global new_vcf

    gz_flag = False
    if is_file_gzipped(filename):
        vcf_file = gzip.open(filename,"r")
        gz_flag = True
    else:
        vcf_file = open(filename,"r") 

    for line in vcf_file:
        header_line = ""
        if not gz_flag:
            fields = line.rstrip().split('\t')
        else:
            header_line = line.decode("utf-8").rstrip()
            fields = header_line.split('\t')
        if fields[0][0] == "#":
            new_vcf += header_line + "\n"
            continue
        
        chr_name = fields[0].split(' ')[0]
        start_pos = int(fields[1]) - 1 # convert to 0-based
        var_id = fields[2]
        ref = fields[3]
        alt = fields[4]
        qual = fields[5]
        filt = fields[6]
        info = fields[7].rstrip().split(';')
        svtype =  info[1][7:]
        tmp_contig = ""
        contig_id = []
        if svtype == "SNV":
            svlen = 1
            tmp_contig = info[2][11:].rstrip().split(',')
        else:
            svlen = abs(int(info[2][6:]))
            tmp_contig = info[3][11:].rstrip().split(',')

        for i in tmp_contig:
            contig_id.append(i.rstrip().split(':')[0])
        
        genotype = fields[9]
 
        yield VCF(chr_name, start_pos, var_id, ref, alt, qual, filt, svtype, svlen, contig_id, genotype)
    return


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-g","--genome", dest="genome_file",
                        help="Path to the reference genome file in fasta format.")
    parser.add_argument("-v","--vcf", dest="vcf_file",
                        help="Path to the VCF file")

    args = parser.parse_args()
    ref_genome = pysam.Fastafile(args.genome_file)

    vars_by_chr_tmp = {}
    vars_by_chr = {}
    dups_by_chr = {}

    for cnt, line in enumerate(parse_vcf(args.vcf_file)):
        if line.svlen < 50:
            continue
        if line.chr_name not in vars_by_chr_tmp:
            vars_by_chr_tmp[line.chr_name] = [line]
        else:
            vars_by_chr_tmp[line.chr_name].append(line)

    #Add only the non duplicates from the same SVtig (contig_id)
    for chr_name, chr_list in vars_by_chr_tmp.items():
        for var1 in chr_list:
            toAdd = True
            for var2 in chr_list:
                if var1 == var2:
                    continue

                if lists_overlap(var1.contig_id, var2.contig_id) and is_overlapping(var1.start_pos, var1.end_pos, var2.start_pos, var2.end_pos) == True:
                    var1.hasDups = True
                    if var1.genotype != var2.genotype:
                        if (var1.genotype == "0|1" or var1.genotype == "1|0") and var2.genotype == "1|1":
                            toAdd = False
                            break
                    if var1.svlen < var2.svlen:
                        toAdd = False
                        break
            if toAdd:
                if var1.chr_name not in vars_by_chr:
                    vars_by_chr[var1.chr_name] = [var1]
                else:
                    vars_by_chr[var1.chr_name].append(var1)
            else:
                if var1.chr_name not in dups_by_chr:
                    dups_by_chr[var1.chr_name] = [var1]
                else:
                    dups_by_chr[var1.chr_name].append(var1)

    contigs = {}
    for k, v in vars_by_chr.items():
        for v2 in v:
            for c in v2.contig_id:
                if c not in contigs:
                    tmp_contig = Contig(k)
                    tmp_contig.vars.append(v2)
                    tmp_contig.genotype = c.rstrip().split('-')[0]
                    contigs[c] = tmp_contig
                else:
                    contigs[c].vars.append(v2)

    ########################################
    for k, v in contigs.items(): 
        for i in v.vars:
            if i.chr_name != v.chr_name:
                break
            if i.start_pos < v.start_pos:
                v.start_pos = i.start_pos
            if i.start_pos + i.svlen > v.end_pos:
                v.end_pos = i.start_pos + i.svlen
    
        v.ref = ref_genome.fetch(v.chr_name, v.start_pos, v.end_pos + 1)
        seq = v.ref[:]

        #Add the variants to the sequence now 
        diff = 0
        for i in v.vars:
            if i.chr_name != v.chr_name:
                break
            offset = i.start_pos - v.start_pos + diff
            sv_size = abs(len(i.ref) - len(i.alt))

            if len(i.ref) > len(i.alt):
                diff -= sv_size
            elif len(i.ref) < len(i.alt):
                diff += sv_size

            subseq = seq[offset:offset + len(i.ref)]
            if subseq.lower() == i.ref.lower():
                tmp_seq = seq[0:offset] + i.alt + seq[offset + len(i.ref):]
                seq = tmp_seq[:]
        
        var_name = v.chr_name + "-" + str(v.start_pos + 1) + "-" + v.genotype
        if v.genotype == "H1":
            genotype = "1|0"
        else:
            genotype = "0|1"
        new_vcf += v.chr_name + "\t" + str(v.start_pos + 1) + "\t" + var_name + "\t" + v.ref + "\t" + seq + "\t.\t.\t" + "ID=" + var_name + ";SVTYPE=COMPLEX;TIG_REGION=" + k + "\t" + "GT\t" + genotype + "\n"
    
    output_file = args.vcf_file+"_merged.vcf"
    print("New VCF is written to", output_file)
    with open(output_file, "w") as text_file:
        text_file.write(new_vcf)


