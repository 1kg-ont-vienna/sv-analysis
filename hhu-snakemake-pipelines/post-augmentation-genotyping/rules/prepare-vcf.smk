configfile: 'config.yaml'

# creating the list of pseudohaplotypes available for all callsets
pseudohaplotypes = {}
for callset in config['callsets'].keys():
    pseudohaplotypes[callset] = {}
    if config['callsets'][callset]['pseudohaplotypes'] == None:
        break
    for line in open(config['callsets'][callset]['pseudohaplotypes'], 'r'):
        line = line.strip()
        name = line.split('/')[-1].split('.')[0]
        pseudohaplotypes[callset][name] = line

# extract HPRC assemblies from AGC
# TODO: Can add a rule for doing this extraction
# if assemblies have been deleted, there is a script in /gpfs/project/projects/medbioinf/users/spani/files/fasta/HPRC/Assemblies_Yr1/ to extract them.

# align HPRC assemblies
rule align_hprc_assemblies:
    input:
        ref='results/{callset}/bubble_calling_and_tagging/{callset}.tagged.gfa',
        assembly=config['path_to_hprc_assemblies']+'/{sample}.{haplotype}.fa',
        path_to_minigraph=config['path_to_minigraph']
    output:
        'results/{callset}/hprc_assembly_mappings/{sample}.{haplotype}.gaf'
    log:
        'results/{callset}/hprc_assembly_mappings/{sample}.{haplotype}.log'
    wildcard_constraints:
        sample='|'.join(config['hprc_samples'])
    resources:
        runtime_hrs=12,
        runtime_min=0,
        mem_total_mb=80000
    threads: 8
    shell:
        '''
        {input.path_to_minigraph} --version > {log}
        {input.path_to_minigraph} --vc -cx lr {input.ref} {input.assembly} -t {threads} > {output} 2>> {log}
        '''

# align pseudohaplotypes
rule align_pseudohaplotypes:
    input:
        ref='results/{callset}/bubble_calling_and_tagging/{callset}.tagged.gfa',
        assembly=lambda wildcards: pseudohaplotypes[wildcards.callset][wildcards.name],
        path_to_minigraph=config['path_to_minigraph']
    output:
        'results/{callset}/pseudo_assembly_mappings/{name}.gaf'
    log:
        'results/{callset}/pseudo_assembly_mappings/{name}.log'
    resources:
        runtime_hrs=12,
        runtime_min=0,
        mem_total_mb=80000
    threads: 8
    shell:
        '''
        {input.path_to_minigraph} --version > {log}
        {input.path_to_minigraph} --vc -cx lr {input.ref} {input.assembly} -t {threads} > {output} 2>> {log}
        '''


# make a list of HPRC assembly GAFs
rule list_hprc_assembly_mappings:
    input:
        expand('results/{{callset}}/hprc_assembly_mappings/{sample}.{haplotype}.gaf', sample=config['hprc_samples'], haplotype=['1', '2'])
    output:
        'results/{callset}/hprc_assembly_mappings/pathlist.txt'
    run:
        f = open(output[0], 'w')
        for name in input:
            print(name, file=f)
        f.close()  


# make a list of pseudo-haplotype GAFs
rule list_pseudohaplotype_assembly_mappings:
    input:
        lambda wildcards: expand('results/{{callset}}/pseudo_assembly_mappings/{name}.gaf', name=[p for p in pseudohaplotypes[wildcards.callset].keys()])
    output:
        'results/{callset}/pseudo_assembly_mappings/pathlist.txt'
    run:
        f = open(output[0], 'w')
        for name in input:
            print(name, file=f)
        f.close()  


# script to make the VCF taking the two lists made above
rule assemblies_to_vcf:
    input:
        ref='results/{callset}/bubble_calling_and_tagging/{callset}.tagged.gfa',
        hprc_path='results/{callset}/hprc_assembly_mappings/pathlist.txt',
        pseudo_path='results/{callset}/pseudo_assembly_mappings/pathlist.txt'
    output:
        'results/{callset}/panel/{callset}.vcf',
        'results/{callset}/panel/{callset}-giggles.vcf',
    params:
        'results/{callset}/panel/'
    log:
        'results/{callset}/panel/{callset}.log'
    conda:
        '../envs/basic.yml'
    resources:
        runtime_hrs=1,
        runtime_min=0,
        mem_total_mb=50*1024
    shell:
        'python scripts/assembly-to-vcf.py -gfa {input.ref} -hprc-list {input.hprc_path} -pseudo-list {input.pseudo_path} -output {params} 2> {log}'


# check correctness of vcfs
rule vcf_correctness:
    input:
        vcf='results/{callset}/panel/{callset}.vcf',
        ref=config['path_to_hprc_mg_ref']
    output:
        out1='results/{callset}/panel/{callset}.check'
    log:
        log1='results/{callset}/panel/{callset}.check.log'
    conda:
        '../envs/basic.yml'
    shell:
        'bcftools norm --check-ref w -f {input.ref} {input.vcf} > /dev/null 2> {log.log1} && touch {output.out1}'


# filter the VCF
rule filter_vcf:
    input:
        vcf1='results/{callset}/panel/{callset}.vcf',
        vcf2='results/{callset}/panel/{callset}-giggles.vcf',
        chk='results/{callset}/panel/{callset}.check'
    output:
        vcf1='results/{callset}/panel/{callset}_filtered.vcf',
        vcf2='results/{callset}/panel/{callset}-giggles_filtered.vcf'
    wildcard_constraints:
        vcf=lambda wildcards: '|'.join(['%s.vcf|%s-giggles.vcf'%(w,w) for w in config['callsets'].keys()])
    conda:
        '../envs/basic.yml'
    shell:
        '''
        bcftools view -f "PASS" -o {output.vcf1} {input.vcf1}
        bcftools view -f "PASS" -o {output.vcf2} {input.vcf2}
        '''

# prepare vcf panel
rule trim_panel:
    input:
        'results/{callset}/panel/{vcf_name}_filtered.vcf'
    output:
        temp('results/{callset}/panel/{vcf_name}_filtered_trimmed.vcf')
    conda:
        '../envs/basic.yml'
    resources:
        mem_total_mb=20000,
        runtime_hrs=0,
        runtime_min=30
    shell:
        'bcftools view --trim-alt-alleles -c1 {input} > {output}'


# annotate giggles vcf with allele decomposition information
rule annotate_panel_giggles_ready:
    input:
        vcf='results/{callset}/panel/{callset}-giggles_filtered_trimmed.vcf',
        gfa='results/{callset}/bubble_calling_and_tagging/{callset}.tagged.gfa'
    output:
        multi='results/{callset}/panel/giggles-ready_multiallelic.vcf',
        multi_tmp=temp('results/{callset}/panel/giggles-tmp.vcf'),
        biallelic='results/{callset}/panel/giggles-ready_biallelic.vcf',
        bi_tmp=temp('results/{callset}/panel/giggles-tmp_biallelic.vcf')
    log:
        'results/{callset}/panel/giggles-ready.annotate.log'
    conda:
        '../envs/basic.yml'
    resources:
        mem_total_mb=20000,
        runtime_hrs=5,
        runtime_min=59
    params:
        outname='results/{callset}/panel/giggles-tmp'
    shell:
        '''
        python3 scripts/annotate_vcf.py -vcf {input.vcf} -gfa {input.gfa} -o {params.outname} &> {log}
        cat {output.multi_tmp} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | \"sort -k1,1 -k2,2n\"}}' > {output.multi}
        cat {output.bi_tmp} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | \"sort -k1,1 -k2,2n\"}}' > {output.biallelic}
        '''

# annotate full vcf with allele decomposition information
rule annotate_panel_full:
    input:
        vcf='results/{callset}/panel/{callset}_filtered_trimmed.vcf',
        gfa='results/{callset}/bubble_calling_and_tagging/{callset}.tagged.gfa'
    output:
        multi='results/{callset}/panel/panel-full_multiallelic.vcf',
        multi_tmp=temp('results/{callset}/panel/panel-tmp.vcf'),
        biallelic='results/{callset}/panel/panel-full_biallelic.vcf',
        bi_tmp=temp('results/{callset}/panel/panel-tmp_biallelic.vcf')
    log:
        'results/{callset}/panel/panel-full.annotate.log'
    conda:
        '../envs/basic.yml'
    resources:
        mem_total_mb=20000,
        runtime_hrs=5,
        runtime_min=59
    params:
        outname='results/{callset}/panel/panel-tmp'
    shell:
        '''
        python3 scripts/annotate_vcf.py -vcf {input.vcf} -gfa {input.gfa} -o {params.outname} &> {log}
        cat {output.multi_tmp} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | \"sort -k1,1 -k2,2n\"}}' > {output.multi}
        cat {output.bi_tmp} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | \"sort -k1,1 -k2,2n\"}}' > {output.biallelic}
        '''