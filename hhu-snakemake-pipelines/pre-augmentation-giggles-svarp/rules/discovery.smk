'''
Issue with 9 sample wrt memory and compute time:
- SVarp calls for HG00139, HG00235, HG00268, HG01699, HG01791, HG01879, HG01915, HG03009, HG03521 were run separately by Arda.
- Reason: Gaftools realign step taking too much memory and time which HILBERT was unable to handle.
- HG00268, HG01791, HG01879, HG01915, HG03009, HG03521 have been run using “-s 10 -d 100”. All other SVarp runs were done using “-s 5 -d 500"
- Samples separately ran by Arda and then added to the snakemake pipeline and touched to further proceed with PAV calling and post processing.

Issue with 4 samples wrt PAV calling:
- Samples are HG02661, HG01600, HG01308, NA19676.
- Very low coverage leading to no SV discovery.
- PAV not able to output empty VCF file and hence these samples were hardcoded out of the svarp call list to allow pipeline to run smoothly.
'''

# run svarp
rule svarp:
    input:
        reads='result/svarp-giggles/data/fasta/{sample}.fasta.gz',
        reads_gzi='result/svarp-giggles/data/fasta/{sample}.fasta.gz.gzi',
        reads_fai='result/svarp-giggles/data/fasta/{sample}.fasta.gz.fai',
        haplotag=config['path_to_haplotags']+'/{sample}/{sample}.tsv',
        alignment='result/svarp-giggles/data/gaf/{sample}.gaf',
        gfa='results/annotating-paths/bubble_calling_and_tagging/chm13-90c.r518_tagged_withseq.gfa'
    output:
        'result/svarp-giggles/svarp/{sample}/{sample}_svtigs_H1.fa',
        'result/svarp-giggles/svarp/{sample}/{sample}_svtigs_H2.fa',
        'result/svarp-giggles/svarp/{sample}/{sample}_svtigs_untagged.fa',
        path_to_svarp=config['path_to_svarp']
    params:
        outdir='result/svarp-giggles/svarp/{sample}',
        path_to_wtdbg2=config['path_to_wtdbg2']
    log:
        stdout='result/svarp-giggles/svarp/log/{sample}.stdout',
        stderr='result/svarp-giggles/svarp/log/{sample}.stderr'
    conda:
        "../envs/svarp.yml"
    resources:
        mem_total_mb=1024*400,
        runtime_hrs=24*5,
        runtime_min=1
    priority: 2
    shell:
        '''
        export PATH=$PATH:{params.path_to_wtdbg2}
        {input.path_to_svarp} -a {input.alignment} -g {input.gfa} --fasta {input.reads} -i {wildcards.sample} --phase {input.haplotag} -o {params.outdir} 2> {log.stderr} 1> {log.stdout}
        '''

# map svtigs with minimap2
rule svtigs_minimap2:
    input:
        ref=config['reference_directory']+'/1KG_ONT_VIENNA_hg38.fa',
        fasta='result/svarp-giggles/svarp/{sample}/{sample}_svtigs_{haplotype}.fa'
    output:
        temp('result/svarp-giggles/svarp/{sample}/svimasm/{sample}_{haplotype}.sam')
    wildcard_constraints:
        haplotype='H1|H2'
    conda:
        '../envs/svarp_processing.yml'
    resources:
        mem_total_mb=30000
    threads: 2
    priority:3
    shell:
        'minimap2 -a -x asm5 --cs -r2k -t {threads} {input.ref} {input.fasta}  > {output}'

# sort aligned SAM file
rule sort_svtig_alignment:
    input:
        'result/svarp-giggles/svarp/{sample}/svimasm/{sample}_{haplotype}.sam'
    output:
        temp('result/svarp-giggles/svarp/{sample}/svimasm/{sample}_{haplotype}.sorted.bam')
    wildcard_constraints:
        haplotype='H1|H2'
    conda:
        '../envs/svarp_processing.yml'
    threads: 2
    priority:3
    shell:
        'samtools sort -m4G -@{threads} -o {output} {input}'

# index sorted BAM
rule index_sorted_svtig_alignment:
    input:
        'result/svarp-giggles/svarp/{sample}/svimasm/{sample}_{haplotype}.sorted.bam'
    output:
        temp('result/svarp-giggles/svarp/{sample}/svimasm/{sample}_{haplotype}.sorted.bam.bai')
    wildcard_constraints:
        haplotype='H1|H2'
    conda:
        '../envs/svarp_processing.yml'
    shell:
        'samtools index {input}'

# run svim-asm to get the vcf
rule run_svim_asm:
    input:
        ref=config['reference_directory']+'/1KG_ONT_VIENNA_hg38.fa',
        bam1='result/svarp-giggles/svarp/{sample}/svimasm/{sample}_H1.sorted.bam',
        bam1_index='result/svarp-giggles/svarp/{sample}/svimasm/{sample}_H1.sorted.bam.bai',
        bam2='result/svarp-giggles/svarp/{sample}/svimasm/{sample}_H2.sorted.bam',
        bam2_index='result/svarp-giggles/svarp/{sample}/svimasm/{sample}_H2.sorted.bam.bai'
    output:
        'result/svarp-giggles/svarp/{sample}/svimasm/variants.vcf',
        'result/svarp-giggles/svarp/{sample}/svimasm/sv-lengths.png'
    params:
        outdir='result/svarp-giggles/svarp/{sample}/svimasm/'
    conda:
        '../envs/svarp_processing.yml'
    priority: 3
    shell:
        'svim-asm diploid --sample {wildcards.sample} {params.outdir} {input.bam1} {input.bam2} {input.ref}'

rule pav_pipeline:
    input:
        'result/svarp-giggles/svarp/{sample}/{sample}_svtigs_H1.fa',
        'result/svarp-giggles/svarp/{sample}/{sample}_svtigs_H2.fa',
        path_to_pav=config['path_to_pav']
    output:
        'result/svarp-giggles/svarp/{sample}/pav_hg38/run.complete',
        'result/svarp-giggles/svarp/{sample}/pav_hg38/pav_svtigs.vcf.gz',
        'result/svarp-giggles/svarp/{sample}/pav_t2t/run.complete',
        'result/svarp-giggles/svarp/{sample}/pav_t2t/pav_svtigs.vcf.gz'
    params:
        dir='result/svarp-giggles/svarp/{sample}/',
        path_to_ref=config['reference_directory']
    threads: 8
    resources:
        runtime_hrs=12,
        mem_total_mb=70000
    shell:
        '''
        cd {params.dir}
        mkdir -p pav_hg38
        mkdir -p pav_t2t
        snakemake -s rules/pav.smk -c {threads} --use-singularity --config sample={wildcards.sample} --config path_to_pav={input.path_to_pav} --config path_to_ref={params.path_to_ref} -F
        rm -r pav_hg38/data
        rm -r pav_hg38/temp
        rm -r pav_hg38/results
        rm -r pav_t2t/data
        rm -r pav_t2t/temp
        rm -r pav_t2t/results
        '''

rule process_pav_output:
    input:
        ref=config['reference_directory']+'/1KG_ONT_VIENNA_{ref}.fa',
        vcf='result/svarp-giggles/svarp/{sample}/pav_{ref}/pav_svtigs.vcf.gz'
    output:
        final='result/svarp-giggles/svarp/{sample}/pav_{ref}/pav_svtigs.vcf.gz_merged.vcf'
    conda:
        '../envs/svarp_processing.yml'
    resources:
        mem_total_mb=5000
    wildcard_constraints:
        ref='t2t|hg38'
    shell:
        'python scripts/svtig_to_single_contig.py -g {input.ref} -v {input.vcf}'

rule rename_pav_output:
    input:
        'result/svarp-giggles/svarp/{sample}/pav_{ref}/pav_svtigs.vcf.gz_merged.vcf'
    output:
        'result/svarp-giggles/svarp/{sample}/pav_{ref}/pav_svtigs_merged.vcf'
    wildcard_constraints:
        ref='t2t|hg38'
    resources:
        mem_total_mb=5000
    shell:
        'mv {input} {output}'