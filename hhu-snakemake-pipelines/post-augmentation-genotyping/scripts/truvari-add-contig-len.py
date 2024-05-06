import sys

genome_index = sys.argv[1]

chromosomes = ['chr' + str(i) for i in range(1,23)]

chrom_to_length = {}
for line in open(genome_index, 'r'):
    fields = line.split()
    if fields[0] in chromosomes:
        chrom_to_length[fields[0]] = fields[1]

counter = 0

no_contigs = True
for line in sys.stdin:
    if line.startswith('##'):
        if line.startswith('##contig=<'):
            no_contigs = False
            # add chromosome length, this is required by truvari..
            f = {s.split('=')[0] : s.split('=')[1] for s in line.strip()[10:-1].split(',') }
            assert 'ID' in f
            chrom = f['ID']
            if not chrom in chromosomes:
                continue
            if not 'length' in f:
                assert chrom in chrom_to_length
                line = "##contig=<ID=" + chrom + ',length=' + chrom_to_length[chrom] + '>'
        print(line.strip())
        continue
    if line.startswith('#'):
        if no_contigs:
            for chrom,length in chrom_to_length.items():
                print("##contig=<ID=" + chrom + ",length=" + length + ">")
        print(line.strip())
        continue
    fields = line.split()
    if fields[0] not in chromosomes:
        continue
    print(line.strip())