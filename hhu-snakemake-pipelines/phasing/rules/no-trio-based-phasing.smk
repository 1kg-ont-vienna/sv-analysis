rule phase_sample_longread:
    input:
        vcf='results/data/nygc-genotypes/no-trios/{sample}/{chr}_filtered.vcf',
        bam=config['path_to_cram']+'/{sample}.hg38.cram',
        ref=config['reference_directory']+'/1KG_ONT_VIENNA_hg38.fa'
    output:
        vcf=temp('results/phased-vcf/longread/{sample}/{chr}.vcf'),
        out='results/phased-vcf/longread/{sample}/{chr}.out'
    conda:
        '../envs/whatshap.yaml'
    resources:
        runtime_hrs=lambda wildcards, attempt: attempt,
        runtime_min=0,
        mem_total_mb=lambda wildcards, attempt: 2048 * attempt,
        mem_per_cpu_mb=lambda wildcards, attempt: 2048 * attempt
    shell:
        '''
        whatshap phase -o {output.vcf} --chromosome {wildcards.chr} --sample {wildcards.sample} -r {input.ref} {input.vcf} {input.bam} 2>&1 | tee {output.out} 
        '''

rule concat_longread:
    input:
        expand('results/phased-vcf/longread/{{sample}}/{chr}.vcf', chr=config['chromosome'])
    output:
        'results/phased-vcf/longread/{sample}.vcf'
    wildcard_constraints:
        chr='[c][h][r][0-9X]{1,2}',
        sample='(?:NA|HG)\d{5}',
        vtype='[a-z]{3,5}'
    conda:
        '../envs/basic.yaml'
    resources:
        runtime_hrs=lambda wildcards, attempt: 3 * attempt,
        runtime_min=0,
        mem_total_mb=lambda wildcards, attempt: 512 * attempt,
        mem_per_cpu_mb=lambda wildcards, attempt: 512 * attempt
    shell:
        '''
        bcftools concat -o {output} {input}
        '''

rule stats_longread:
    input:
        'results/phased-vcf/longread/{sample}.vcf'
    output:
        'results/phased-vcf/longread/{sample}.stats.tsv'
    conda:
        '../envs/whatshap.yaml'
    resources:
        runtime_hrs=0,
        runtime_min=30,
        mem_total_mb=lambda wildcards, attempt: 2048 * attempt,
    shell:
        'whatshap stats --tsv={output} {input}' 