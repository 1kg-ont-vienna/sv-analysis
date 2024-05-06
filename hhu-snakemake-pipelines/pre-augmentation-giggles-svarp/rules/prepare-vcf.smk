chromosomes = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM']

# Global Wildcard Constraints
wildcard_constraints:
    sample='(?:NA|HG)(?:\d{5}|\d{3})',
    chr='chr[0-9A-Z]',
    haplotype='(1)|(2)'

################################################
##### Aligning the Assemblies to the Graph #####
################################################

rule align_assemblies_to_chm13_minigraph:
    input:
        ref=config['path_to_hprc_mg'],
        assembly=config['path_to_hprc_assemblies']+'/{sample}.{haplotype}.fa',
        path_to_minigraph=config['path_to_minigraph']
    output:
        'results/annotating-paths/assembly_mappings/{sample}.{haplotype}.gaf'
    log:
        'results/annotating-paths/log/assembly_mappings/{sample}.{haplotype}.log'
    resources:
        runtime_hrs=12,
        runtime_min=0,
        mem_total_mb=lambda wildcards, attempt: 80000 * attempt
    threads: 24
    shell:
        '{input.path_to_minigraph} --vc -cx lr {input.ref} {input.assembly} -t {threads} > {output}'


rule make_assembly_mapping_list:
    input:
        expand('results/annotating-paths/assembly_mappings/{sample}.{haplotype}.gaf', sample=config['samples'], haplotype=config['haplotypes'])
    output:
        'results/annotating-paths/assembly_mappings/pathlist.txt'
    resources:
        runtime_hrs=0,
        runtime_min=10,
        mem_total_mb=lambda wildcards, attempt: 2048 * attempt,
        mem_per_cpu_mb=lambda wildcards, attempt: 2048 * attempt
    run:
        f = open(output[0], 'w')
        for name in input:
            print(name, file=f)
        f.close()    


rule get_alignment_stats:
    input:
        assembly_list='results/annotating-paths/assembly_mappings/pathlist.txt',
        script='scripts/alignment_stats.py'
    output:
        pics=expand('results/annotating-paths/assembly_mappings/stats/{sample}_{haplotype}_{graph}.png', sample=config['samples'], haplotype=config['haplotypes'], graph=['dv', 'mapq', 'num_align', 'num_res_match', 'p_align']),
        txt='results/annotating-paths/assembly_mappings/stats/text_stats.txt'
    params:
        'results/annotating-paths/assembly_mappings/stats/'
    conda:
        '../envs/basic.yml'
    resources:
        runtime_hrs=0,
        runtime_min=30,
        mem_total_mb=lambda wildcards, attempt: 20*1024 * attempt
    shell:
        'python {input.script} --assembly-list {input.assembly_list} --output {params} > {output.txt}'


############################################
##### Pre-processing GFA and GAF files #####
############################################

# conda environment should have networkx and whatshap. In this version, the tagged GFA does not have sequence info.
rule call_rGFA_bubbles:
    input:
        ref=config['path_to_hprc_mg']
    output:
        expand('results/annotating-paths/bubble_calling_and_tagging/chm13-90c.r518-{chr}.gfa', chr=chromosomes),
        expand('results/annotating-paths/bubble_calling_and_tagging/chm13-90c.r518-{chr}.csv', chr=chromosomes)
    params:
        out_dir='results/annotating-paths/bubble_calling_and_tagging'
    conda:
        '../envs/basic.yml'
    resources:
        runtime_hrs=0,
        runtime_min=30,
        mem_total_mb=lambda wildcards, attempt: 10*1024 * attempt
    shell:
        '''
        set +u
        source ~/.bashrc
        conda activate gaftools-env
        set -u
        gaftools order_gfa --with-sequence --outdir {params.out_dir} {input.ref}
        '''

#Concat the chromsome-wise tagged GFA produced
rule concat_tagged_GFA:
    input:
        expand('results/annotating-paths/bubble_calling_and_tagging/chm13-90c.r518-{chr}.gfa', chr=chromosomes)
    output:
        'results/annotating-paths/bubble_calling_and_tagging/chm13-90c.r518_tagged_withseq.gfa'
    shell:
        'cat {input} > {output}'

# Sorting GAF generated from Assembly Alignment
rule sort_GAF:
    input:
        gfa='results/annotating-paths/bubble_calling_and_tagging/chm13-90c.r518_tagged_withseq.gfa',
        gaf='results/annotating-paths/assembly_mappings/{sample}.{haplotype}.gaf',
    output:
        'results/annotating-paths/assembly_mappings/{sample}.{haplotype}.sorted.gaf'
    log:
        'results/annotating-paths/log/assembly_mappings/{sample}.{haplotype}.sorted.log'
    conda:
        '../envs/basic.yml'
    resources:
        runtime_hrs=1,
        runtime_min=0,
        mem_total_mb=lambda wildcards, attempt: 2048 * attempt
    shell:
        '''
        set +u
        source ~/.bashrc
        conda activate gaftools-env
        set -u
        gaftools sort --bgzip --outgaf {output} {input.gaf} {input.gfa} 2> {log}
        '''
    

##################################
##### Producing the VCF file #####
##################################

rule assemblies_to_vcf:
    input:
        gfa='results/annotating-paths/bubble_calling_and_tagging/chm13-90c.r518_tagged_withseq.gfa',
        assembly_list='results/annotating-paths/assembly_mappings/pathlist.txt',
        script='scripts/assembly-to-vcf.py'
    output:
        'results/annotating-paths/vcf/chm13-90c.r518.vcf.gz'
    log:
        'results/annotating-paths/log/vcf/chm13-90c.r518.log'
    conda:
        '../envs/basic.yml'
    resources:
        runtime_hrs=1,
        runtime_min=0,
        mem_total_mb=lambda wildcards, attempt: 10*1024 * attempt
    shell:
        '''
        python {input.script} --gfa {input.gfa} --assembly-list {input.assembly_list} --bgzip --output {output} 2> {log}
        bcftools index {output}
        '''

# Filter out the NS80 tagged records
rule filtered_vcf:
    input:
        vcf='results/annotating-paths/vcf/chm13-90c.r518.vcf.gz',
        check='results/annotating-paths/vcf/chm13-90c.r518.check'
    output:
        'results/annotating-paths/vcf/chm13-90c.r518_filtered.vcf.gz'
    conda:
        '../envs/basic.yml'
    resources:
        runtime_hrs=1,
        runtime_min=0,
        mem_total_mb=lambda wildcards, attempt: 2*1024 * attempt
    shell:
        '''
        bcftools view -f "PASS" -Oz -o {output} {input}
        bcftools index {output}
        '''

# check correctness of vcfs
rule vcf_correctness:
    input:
        vcf='results/annotating-paths/vcf/chm13-90c.r518.vcf.gz',
        ref=config['path_to_hprc_mg_ref']
    output:
        out1='results/annotating-paths/vcf/chm13-90c.r518.check'
    log:
        log1='results/annotating-paths/vcf/chm13-90c.r518.check.log'
    conda:
        '../envs/basic.yml'
    shell:
        'bcftools norm --check-ref w -f {input.ref} {input.vcf} > /dev/null 2> {log.log1} && touch {output.out1}'

# prepare vcf panel
rule prepare_vcf:
    input:
        'results/annotating-paths/vcf/chm13-90c.r518_filtered.vcf.gz'
    output:
        'result/svarp-giggles/data/vcf/panel.vcf'
    conda:
        "../envs/basic.yml"
    resources:
        mem_total_mb=20000,
        runtime_hrs=0,
        runtime_min=30
    shell:
        "bcftools view --trim-alt-alleles -c1 {input} > {output}"