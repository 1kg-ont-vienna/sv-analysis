rule no_trio_comparison:
    input:
        longread='results/phased-vcf/longread/{sample}.vcf',
        nygc='results/data/nygc-phased/no-trios/{sample}/filtered.vcf'
    output:
        pairwise='results/no-trio-comparison/{sample}/pairwise.tsv',
        multiway='results/no-trio-comparison/{sample}/multiway.tsv'
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
        whatshap compare --sample {wildcards.sample} --only-snvs --names longread,nygc --tsv-pairwise {output.pairwise} --tsv-multiway {output.multiway} {input.longread} {input.nygc}
        '''
