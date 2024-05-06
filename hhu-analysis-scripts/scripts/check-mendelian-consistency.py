import argparse
import gzip
import copy

def check_mendelian_consistency(child_gt, parent1_gt, parent2_gt):
    child = child_gt.get_alleles()
    parent1 = parent1_gt.get_alleles()
    parent2 = parent2_gt.get_alleles()
    if child[0] in parent1 and child[1] in parent2:
        return True
    if child[0] in parent2 and child[1] in parent1:
        return True
    return False

class Genotype:
    def __init__(self, alleles, is_phased):
        self._alleles = alleles
        self._phased = is_phased

    def __str__(self):
        if self._phased:
            return '|'.join([str(i) for i in self._alleles])
        else:
            return '/'.join([str(i) for i in self._alleles])

    def get_alleles(self):
        return self._alleles
    
    def get_index(self):
        return sum(self._alleles)

    def get_ploidy(self):
        return len(self._alleles)

    def is_phased(self):
        return self._phased

    def __eq__(self, other):
        return sorted(self._alleles) == sorted(other._alleles)

    def is_hom_ref(self):
        return all([int(a) == 0 for a in self._alleles])

    def is_none(self):
        return self._alleles == []

def genotype_from_string(gt_string):
    is_phased = False
    alleles = []
    if '.' in gt_string:
        # untyped
        return Genotype(alleles, is_phased)

    if '|' in gt_string:
        is_phased = True
        alleles = [int(allele) for allele in gt_string.split('|')]
    elif '/' in gt_string:
        is_phased = False
        alleles = [int(allele) for allele in gt_string.split('/')]
    else:
        assert False
    return Genotype(alleles, is_phased)


def run_statistics(vcf_file, father, mother, child, out, var_type, bub2type, allele2bub):
    bed_writer = open(out+'/'+child+'-'+var_type+'.bed', 'w')
    table_writer = open(out+'/'+child+'-'+var_type+'.tsv', 'w')
    print('father_gt\tmother_gt\tchild_gt\tcount\tconsistent', file=table_writer)
    sample_to_index = {}
    count = {}
    for x in [0, 1, 2]:
        for y in [0, 1, 2]:
            for z in [0, 1, 2]:
                count[(x, y, z)] = 0
    
    total_consistent_trios = 0
    alt_transmitted = 0
    all_het = 0
    all_abs = 0
    all_present = 0
    consistent_trios = 0
    inconsistent_trios = 0
    untyped_trios = 0
    total_untrivial_trios = 0

    total_consistent_trios_vtype = {'COMPLEX': 0, 'DEL': 0, 'INS': 0, 'SNV': 0}
    alt_transmitted_vtype = copy.deepcopy(total_consistent_trios_vtype)
    all_het_vtype = copy.deepcopy(total_consistent_trios_vtype)
    all_abs_vtype = copy.deepcopy(total_consistent_trios_vtype)
    all_present_vtype = copy.deepcopy(total_consistent_trios_vtype)
    consistent_trios_vtype = copy.deepcopy(total_consistent_trios_vtype)
    inconsistent_trios_vtype = copy.deepcopy(total_consistent_trios_vtype)
    untyped_trios_vtype = copy.deepcopy(total_consistent_trios_vtype)
    total_untrivial_trios_vtype = copy.deepcopy(total_consistent_trios_vtype)


    for record in gzip.open(vcf_file, 'rt'):
        if record.startswith('##'):
            continue
        if record.startswith('#'):
            header = record.strip().split()
            for i,f in enumerate(header):
                sample_to_index[f] = i
            continue
        assert header is not None
        

        fields = record.split()
        genotype_index = fields[8].split(':').index('GT')
        assert len(fields[4].split(',')) == 1
        info_fields = { i.split('=')[0] : i.split('=')[1].strip() for i in fields[7].split(';') if '=' in i} 
        assert 'ID' in info_fields
        assert len(info_fields['ID'].split(',')) == 1
        variant_id = info_fields['ID']
        variant_chrom = fields[0]
        variant_start = fields[1]
        variant_end = str(int(fields[1]) + len(fields[3]))
        variant_type = variant_id.split('-')[2]
        
        if var_type == 'all':
            pass
        elif var_type == 'biallelic':
            if bub2type[allele2bub[variant_id]]:
                continue
        elif var_type == 'multiallelic':
            if not bub2type[allele2bub[variant_id]]:
                continue
        
        if 'X' in fields[0] or 'Y' in fields[0]:
            continue
        field_child = fields[sample_to_index[child]].split(':')
        field_father = fields[sample_to_index[father]].split(':')
        field_mother = fields[sample_to_index[mother]].split(':')
        gt_child = genotype_from_string(field_child[genotype_index])
        gt_father = genotype_from_string(field_father[genotype_index])
        gt_mother = genotype_from_string(field_mother[genotype_index])
        if any([g.is_none() for g in [gt_child, gt_father, gt_mother]]):
            untyped_trios += 1
            untyped_trios_vtype[variant_type] += 1
        elif all([ g == gt_child for g in [gt_father, gt_mother]]):
            # all genotypes same, automatically consistent
            if gt_child == Genotype([0,0], False):
                all_abs += 1
                all_abs_vtype[variant_type] += 1
            elif gt_child == Genotype([0,1], False):
                all_het += 1
                all_het_vtype[variant_type] += 1
            elif gt_child == Genotype([1,1], False):
                all_present += 1
                all_present_vtype[variant_type] += 1
            else:
                assert(False)
            total_consistent_trios += 1
            total_consistent_trios_vtype[variant_type] += 1
            count[(gt_father.get_index(), gt_mother.get_index(), gt_child.get_index())] += 1
        else:
            count[(gt_father.get_index(), gt_mother.get_index(), gt_child.get_index())] += 1
            total_untrivial_trios += 1
            total_untrivial_trios_vtype[variant_type] += 1
            consistent = check_mendelian_consistency(gt_child, gt_father, gt_mother)
            if consistent:
                consistent_trios += 1
                consistent_trios_vtype[variant_type] += 1
                total_consistent_trios += 1
                total_consistent_trios_vtype[variant_type] += 1
                # check how often alt allele was transmitted
                alt_transmitted += sum(a!=0 for a in gt_child.get_alleles())
                alt_transmitted_vtype[variant_type] += sum(a!=0 for a in gt_child.get_alleles())
            else:
                inconsistent_trios += 1
                inconsistent_trios_vtype[variant_type] += 1
                print(variant_chrom+'\t'+variant_start+'\t'+variant_end+'\t'+variant_id,file=bed_writer)
    
    index_to_allele_map = {0: [0,0], 1: [0,1], 2: [1,1]}
    for x in [0, 1, 2]:
        for y in [0, 1, 2]:
            for z in [0, 1, 2]:
                is_consistent = check_mendelian_consistency(
                    Genotype(index_to_allele_map[z], False),
                    Genotype(index_to_allele_map[x], False),
                    Genotype(index_to_allele_map[y], False),
                    )
                print('%d\t%d\t%d\t%d\t%s'%(x,y,z,count[x,y,z], str(is_consistent)), file=table_writer)
    
    table_writer.close()
    bed_writer.close()
    
    print("\n#Variant Type:", var_type)    
    print("\n#Statistics:")
    print("\tall\tINS\tDEL\tOTHERS")
    print("Total consistent trios", total_consistent_trios, '\t'.join([str(total_consistent_trios_vtype[x]) for x in ['INS', 'DEL', 'COMPLEX']]), sep='\t')
    print("Total untrivial consistent trios", total_untrivial_trios, '\t'.join([str(total_untrivial_trios_vtype[x]) for x in ['INS', 'DEL', 'COMPLEX']]), sep='\t')
    print("Total inconsistent trios", inconsistent_trios, '\t'.join([str(inconsistent_trios_vtype[x]) for x in ['INS', 'DEL', 'COMPLEX']]), sep='\t')
    print("Untyped trios", untyped_trios, '\t'.join([str(untyped_trios_vtype[x]) for x in ['INS', 'DEL', 'COMPLEX']]), sep='\t')
    print("Alternative allele transmitted", alt_transmitted, '\t'.join([str(alt_transmitted_vtype[x]) for x in ['INS', 'DEL', 'COMPLEX']]), sep='\t')
    print("All ref trios", all_abs, '\t'.join([str(all_abs_vtype[x]) for x in ['INS', 'DEL', 'COMPLEX']]), sep='\t')
    print("All het trios", all_het, '\t'.join([str(all_het_vtype[x]) for x in ['INS', 'DEL', 'COMPLEX']]), sep='\t')
    print("All hom trios", all_present, '\t'.join([str(all_present_vtype[x]) for x in ['INS', 'DEL', 'COMPLEX']]), sep='\t')

def plot_statistics(child, out, var_type):
    '''
    plots the confusion matrix as a table
    '''
    column_writer = open(out+'/'+child+'-'+var_type+'.out', 'w')
    table_reader = open(out+'/'+child+'-'+var_type+'.tsv', 'r')
    table = {}
    for line in table_reader:
        if line[0] == 'f':
            continue
        f_gt, m_gt, c_gt, count, consistent = line.rstrip().split('\t')
        c = True if consistent == 'True' else False
        table[(int(f_gt), int(m_gt), int(c_gt))] = [int(count), c]
    table_reader.close()
    
    child = [0,1,2,0,1,2,0,1,2]
    columns = []
    for f in [0,1,2]:
        tmp = []
        for m in [0,1,2]:
            for c in [0,1,2]:
                tmp.append(table[f,m,c][0])
        columns.append(tmp)
    for i in range(9):
        line = "%d\t%d\t%d\t%d\t%d\t%d"%(child[i], columns[0][i], child[i], columns[1][i], child[i], columns[2][i])
        print(line, file=column_writer)
    column_writer.close()

def read_map(file):
    bub2type={}
    allele2bub={}
    with open(file, 'r') as f:
        for line in f:
            if line[0] == '#':
                continue
            bub_id, allele_ids = line.rstrip().split('\t')
            allele_ids = allele_ids.split(',')
            if len(allele_ids) > 1:
                bub2type[bub_id] = True
            else:
                bub2type[bub_id] = False
            for a in allele_ids:
                allele2bub[a] = bub_id
    return bub2type, allele2bub

parser = argparse.ArgumentParser(prog='mendelian-consistency.py', description=__doc__)
parser.add_argument('-vcf', metavar='VCF', required=True, help='Multisample VCF-file with genotypes.')
parser.add_argument('-father', metavar='FATHER', required=True, help='Father sample')
parser.add_argument('-mother', metavar='MOTHER', required=True, help='Mother sample')
parser.add_argument('-child', metavar='CHILD', required=True, help='Child sample')
parser.add_argument('-map', metavar='MAP', required=True, help='bubble id to allele id map')
parser.add_argument('-out', metavar='OUT', required=True, help='output directory')

args = parser.parse_args()
bub2type, allele2bub =  read_map(args.map)

print("\n#Family:")
print("Father:\t", args.father)
print("Mother:\t", args.mother)
print("Child:\t", args.child)

for var_type in ['all', 'multiallelic', 'biallelic']:
    run_statistics(args.vcf, args.father, args.mother, args.child, args.out, var_type, bub2type, allele2bub)
    plot_statistics(args.child, args.out, var_type)
