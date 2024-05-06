include: './get-sample-list.smk'

wildcard_constraints:
    sample='|'.join(samples),
    pop='|'.join(['AFR', 'AMR', 'EAS', 'EUR', 'SAS']),
    sv_type='|'.join(['COMPLEX', 'DEL', 'INS']),
    ranges='|'.join(['0-1kbp', '1kbp-10kbp', '10kbp-100kbp', '100kbp-1Mbp']),
    regions='|'.join(['biallelic', 'multiallelic']),
    vartype='|'.join(['all', 'sv', 'indels', 'large-deletion', 'large-insertion', 'large-complex'])

# run genotyping
rule giggles:
    input:
        reads='results/data/fasta/{sample}.fasta.gz',
        fai='results/data/fasta/{sample}.fasta.gz.fai',
        gzi='results/data/fasta/{sample}.fasta.gz.gzi',
        haplotag=config['path_to_haplotags']+'/{sample}/{sample}.tsv',
        alignment='results/{callset}/ont-alignments/{sample}.sorted.gaf.gz',
        alignment_index='results/{callset}/ont-alignments/{sample}.sorted.gaf.gz.gai',
        gfa='results/{callset}/bubble_calling_and_tagging/{callset}.tagged.gfa',
        vcf='results/{callset}/panel/giggles-ready_multiallelic.vcf.gz'
    output:
        vcf=temp('results/{callset}/genotypes/{sample}-multiallelic.vcf')
    log:
        stderr='results/{callset}/genotypes/{sample}-multiallelic.stderr'
    resources:
        mem_total_mb=20000,
        runtime_hrs=71,
        runtime_min=59
    shell:
        """
        set +u
        source ~/.bashrc
        conda activate giggles-env
        set -u
        giggles --version > {log}
        giggles genotype --read-fasta {input.reads} --sample {wildcards.sample} --realign-mode edit -o {output.vcf} --rgfa {input.gfa} --haplotag-tsv {input.haplotag} {input.vcf} {input.alignment} 2>> {log.stderr}
        """


# get the biallelic and multiallelic VCFs for the genotypes
rule convert_sample_vcf_biallelic:
    input:
        sample_vcf='results/{callset}/genotypes/{sample}-multiallelic.vcf',
        biallelic_vcf='results/{callset}/panel/giggles-ready_biallelic.vcf.gz'
    output:
        temp('results/{callset}/genotypes/{sample}-biallelic.vcf')
    conda:
        '../envs/basic.yml'
    resources:
        mem_total_mb=3000
    shell:
        "cat {input.sample_vcf} | python scripts/convert-to-biallelic.py {input.biallelic_vcf} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | \"sort -k1,1 -k2,2n \"}}' > {output}"


# merge multiallelic records to give multisample vcf
rule merge_vcf_to_multisample:
    input:
        expand('results/{{callset}}/genotypes/{sample}-multiallelic.vcf.gz', sample=samples)
    output:
        temp('results/{callset}/genotypes/multisample-multiallelic.vcf')
    resources:
        mem_total_mb=180000,
        runtime_hrs=71,
    conda:
        '../envs/basic.yml'
    shell:
        '''
        bcftools merge --no-version -o {output} {input}
        '''


# get the biallelic multisample VCF
rule convert_multisample_vcf_biallelic:
    input:
        vcf='results/{callset}/genotypes/multisample-multiallelic.vcf',
        biallelic_vcf='results/{callset}/panel/giggles-ready_biallelic.vcf.gz'
    output:
        temp('results/{callset}/genotypes/multisample-biallelic.vcf')
    conda:
        '../envs/basic.yml'
    resources:
        mem_total_mb=3000
    shell:
        "cat {input.vcf} | python scripts/convert-to-biallelic.py {input.biallelic_vcf} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | \"sort -k1,1 -k2,2n \"}}' > {output}"