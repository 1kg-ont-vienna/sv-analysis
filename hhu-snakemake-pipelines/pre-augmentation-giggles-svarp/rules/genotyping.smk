include: './get-sample-list.smk'

giggles_memory_req={'HG02013': 50000, 'HG00268': 30000, 'NA21113': 30000}

# add allele decomposition info to the panel VCF and get the biallelic and multiallelic panel VCF.
rule annotate_panel:
    input:
        vcf='result/svarp-giggles/data/vcf/panel.vcf.gz',
        gfa='results/annotating-paths/bubble_calling_and_tagging/chm13-90c.r518_tagged_withseq.gfa'
    output:
        vcf_temp=temp('result/svarp-giggles/data/vcf/panel-unzipped-tmp.vcf'),
        multi='result/svarp-giggles/data/vcf/panel-multiallelic.vcf',
        multi_tmp=temp('result/svarp-giggles/data/vcf/panel-tmp.vcf'),
        biallelic='result/svarp-giggles/data/vcf/panel-biallelic.vcf',
        bi_tmp=temp('result/svarp-giggles/data/vcf/panel-tmp_biallelic.vcf')
    log:
        "result/svarp-giggles/data/vcf/panel-annotate.log"
    conda:
        '../envs/basic.yml'
    resources:
        mem_total_mb=20000,
        runtime_hrs=5,
        runtime_min=59
    params:
        outname='result/svarp-giggles/data/vcf/panel-tmp'
    shell:
        """
        gzip -d -c {input.vcf} > {output.vcf_temp}
        python3 scripts/annotate_vcf.py -vcf {output.vcf_temp} -gfa {input.gfa} -o {params.outname} &> {log}
        cat {output.multi_tmp} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | \"sort -k1,1 -k2,2n\"}}' > {output.multi}
        cat {output.bi_tmp} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | \"sort -k1,1 -k2,2n\"}}' > {output.biallelic}
        """

# run genotyping
rule giggles:
    input:
        reads='result/svarp-giggles/data/fasta/{sample}.fasta.gz',
        reads_gzi='result/svarp-giggles/data/fasta/{sample}.fasta.gz.gzi',
        reads_fai='result/svarp-giggles/data/fasta/{sample}.fasta.gz.fai',
        haplotag=config['path_to_haplotags']+'/{sample}/{sample}.tsv',
        alignment='result/svarp-giggles/data/gaf/{sample}.sorted.gaf.gz',
        alignment_index='result/svarp-giggles/data/gaf/{sample}.sorted.gaf.gz.gai',
        gfa='results/annotating-paths/bubble_calling_and_tagging/chm13-90c.r518_tagged_withseq.gfa',
        vcf='result/svarp-giggles/data/vcf/panel-multiallelic.vcf'
    output:
        vcf=temp('result/svarp-giggles/genotypes/{sample}-multiallelic.vcf')
    log:
        stderr="result/svarp-giggles/genotypes/{sample}.stderr",
        stdout="result/svarp-giggles/genotypes/{sample}.stdout"
    resources:
        mem_total_mb=lambda wildcards: 20000 if wildcards.sample not in ['HG02013', 'HG00268', 'NA21113'] else giggles_memory_req[wildcards.sample],
        runtime_hrs=71,
        runtime_min=59
    priority: 3
    shell:
        """
        set +u
        source ~/.bashrc
        conda activate giggles-env
        set -u
        giggles genotype --read-fasta {input.reads} --sample {wildcards.sample} --realign-mode edit -o {output.vcf} --rgfa {input.gfa} --haplotag-tsv {input.haplotag} {input.vcf} {input.alignment} > {log.stdout} 2> {log.stderr}
        """    

# get the biallelic and multiallelic VCFs for the genotypes
rule convert_sample_vcf_biallelic:
    input:
        sample_vcf='result/svarp-giggles/genotypes/{sample}-multiallelic.vcf',
        biallelic_vcf='result/svarp-giggles/data/vcf/panel-biallelic.vcf.gz'
    output:
        temp('result/svarp-giggles/genotypes/{sample}-biallelic.vcf')
    conda:
        '../envs/basic.yml'
    resources:
        mem_total_mb=3000
    shell:
        "cat {input.sample_vcf} | python scripts/convert-to-biallelic.py {input.biallelic_vcf} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | \"sort -k1,1 -k2,2n \"}}' > {output}"


# merge multiallelic records to give multisample vcf
rule merge_vcf_to_multisample:
    input:
        expand('result/svarp-giggles/genotypes/{sample}-multiallelic.vcf.gz', sample=samples)
    output:
        temp('result/svarp-giggles/genotypes/multisample-multiallelic.vcf')
    conda:
        '../envs/basic.yml'
    resources:
        mem_total_mb=180000,
        runtime_hrs=71,
    shell:
        'bcftools merge --no-version -o {output} {input}'


# get the biallelic multisample VCF
rule convert_multisample_vcf_biallelic:
    input:
        sample_vcf='result/svarp-giggles/genotypes/multisample-multiallelic.vcf',
        biallelic_vcf='result/svarp-giggles/data/vcf/panel-biallelic.vcf.gz'
    output:
        temp('result/svarp-giggles/genotypes/multisample-biallelic.vcf')
    conda:
        '../envs/basic.yml'
    resources:
        mem_total_mb=3000
    shell:
        "cat {input.sample_vcf} | python scripts/convert-to-biallelic.py {input.biallelic_vcf} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | \"sort -k1,1 -k2,2n \"}}' > {output}"