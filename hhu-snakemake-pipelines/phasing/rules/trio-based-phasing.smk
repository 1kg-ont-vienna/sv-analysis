rule phase_trio:
    input:
        vcf='results/data/nygc-genotypes/trios/{family}/{chr}_filtered.vcf',
        ref=config['reference_directory']+'/1KG_ONT_VIENNA_hg38.fa',
        ped='resources/pedigree.ped'
    output:
        vcf=temp('results/phased-vcf/trio_phase/{family}/{chr}.vcf'),
        out='results/phased-vcf/trio_phase/{family}/{chr}.out'
    conda:
        '../envs/whatshap.yaml'
    resources:
        runtime_hrs=lambda wildcards, attempt: 2,
        mem_total_mb=lambda wildcards, attempt: 2048 * attempt,
        mem_per_cpu_mb=lambda wildcards, attempt: 2048 * attempt
    shell:
        '''
        whatshap phase -o {output.vcf} --chromosome {wildcards.chr} -r {input.ref} --ped {input.ped} {input.vcf} 2>&1 | tee {output.out} 
        '''

rule concat_trio_phase:
    input:
        expand('results/phased-vcf/trio_phase/{{family}}/{chr}.vcf', chr=config['chromosome'])
    output:
        'results/phased-vcf/trio_phase/{family}.vcf'
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

rule stats_trio_phase:
    input:
        'results/phased-vcf/trio_phase/{family}.vcf'
    output:
        'results/phased-vcf/trio_phase/{family}.stats.tsv'
    conda:
        '../envs/whatshap.yaml'
    resources:
        runtime_hrs=0,
        runtime_min=30,
        mem_total_mb=lambda wildcards, attempt: 2048 * attempt,
    shell:
        'whatshap stats --tsv={output} {input}'    


rule longread_trio_phase:
    input:
        vcf='results/data/nygc-genotypes/trios/{family}/{chr}_filtered.vcf',
        ref=config['reference_directory']+'/1KG_ONT_VIENNA_hg38.fa',
        bam= lambda wildcards: expand(config['path_to_cram']+'/{sample}.hg38.cram', sample=wildcards.family.split('_')),
        ped='resources/pedigree.ped'
    output:
        vcf=temp('results/phased-vcf/longread_trio_phase/{family}/{chr}.vcf'),
        out='results/phased-vcf/longread_trio_phase/{family}/{chr}.out'
    conda:
        '../envs/whatshap.yaml'
    resources:
        runtime_hrs=lambda wildcards, attempt: attempt,
        mem_total_mb=lambda wildcards, attempt: 2048 * attempt,
        mem_per_cpu_mb=lambda wildcards, attempt: 2048 * attempt
    shell:
        '''
        whatshap phase -o {output.vcf} --chromosome {wildcards.chr} -r {input.ref} --ped {input.ped} {input.vcf} {input.bam} 2>&1 | tee {output.out} 
        '''

rule concat_longread_trio:
    input:
        expand('results/phased-vcf/longread_trio_phase/{{family}}/{chr}.vcf', chr=config['chromosome'])
    output:
        'results/phased-vcf/longread_trio_phase/{family}.vcf'
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

rule stats_longread_trio:
    input:
        'results/phased-vcf/longread_trio_phase/{family}.vcf'
    output:
        'results/phased-vcf/longread_trio_phase/{family}.stats.tsv'
    conda:
        '../envs/whatshap.yaml'
    resources:
        runtime_hrs=0,
        runtime_min=30,
        mem_total_mb=lambda wildcards, attempt: 2048 * attempt,
    shell:
        'whatshap stats --tsv={output} {input}' 