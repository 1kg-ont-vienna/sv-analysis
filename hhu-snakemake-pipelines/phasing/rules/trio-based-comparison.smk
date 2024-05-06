sample2family={}
for family in config['families']:
    for sample in family.split('_'):
        sample2family[sample] = family

rule trio_comparison:
    input:
        trio=lambda wildcards: 'results/phased-vcf/trio_phase/'+sample2family[wildcards.sample]+'.vcf',
        longread_trio=lambda wildcards: 'results/phased-vcf/longread_trio_phase/'+sample2family[wildcards.sample]+'.vcf',
        longread='results/phased-vcf/longread/{sample}.vcf',
        nygc=lambda wildcards: 'results/data/nygc-phased/trios/'+sample2family[wildcards.sample]+'/filtered.vcf'
    output:
        comp12='results/trio-comparison/{sample}/trio_longread.txt',
        comp13='results/trio-comparison/{sample}/trio_triolongread.txt',
        comp1stat='results/trio-comparison/{sample}/trio_stat.txt',
        pairwise='results/trio-comparison/{sample}/pairwise.tsv',
        multiway='results/trio-comparison/{sample}/multiway.tsv'
    wildcard_constraints:
        sample='(?:NA|HG)\d{5}',
        vtype='[a-z]{3,5}'
    conda:
        '../envs/whatshap.yaml'
    resources:
        runtime_hrs=lambda wildcards, attempt: 3 * attempt,
        mem_total_mb=lambda wildcards, attempt: 10240 * attempt,
        mem_per_cpu_mb=lambda wildcards, attempt: 10240 * attempt
    shell:
        '''
        whatshap compare --sample {wildcards.sample} --names trio,longread {input.trio} {input.longread} > {output.comp12}
        whatshap compare --sample {wildcards.sample} --names trio,trio-longread {input.trio} {input.longread_trio} > {output.comp13}
        whatshap compare --sample {wildcards.sample} --names trio,nygc {input.trio} {input.nygc} > {output.comp1stat}
        whatshap compare --sample {wildcards.sample} --names trio,longread,trio-longread,nygc --tsv-pairwise {output.pairwise} --tsv-multiway {output.multiway} {input.trio} {input.longread} {input.longread_trio} {input.nygc}
        '''
