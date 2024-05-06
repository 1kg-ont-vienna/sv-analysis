import sys
import pandas as pd
import gzip

if __name__ == "__main__":
    filename = sys.argv[1]
    inp_vcf = sys.argv[2]
    outname = sys.argv[3]

    # creating the final table with the filters

    df = pd.read_csv(filename, sep='\t', header=0)
    passed_var_ids = df[df['confidence_level'] > 0]['VARIANT_ID'].to_numpy()


    # subsetting vcf
    out_vcf = open(outname+'filtered-final.vcf', 'w')
    out_vcf_lowq = open(outname+'lowq-final.vcf', 'w')
    vcf_reader = gzip.open(inp_vcf, 'r')
    count = 0
    for line in vcf_reader:
        line = line.decode('utf-8')
        if line[0] == '#':
            out_vcf.write(line)
            out_vcf_lowq.write(line)
            continue
        id = line.split('\t')[2]
        if count == passed_var_ids.__len__():
            break
        if id == passed_var_ids[count]:
            count += 1
            out_vcf.write(line)
        else:
            out_vcf_lowq.write(line)
