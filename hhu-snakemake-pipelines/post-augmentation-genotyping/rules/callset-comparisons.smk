# Currently not using any of this results for analysis. Commenting the entire code.

"""
include: './get-sample-list.smk'

max_af = 1
min_af = 0
chromosomes='chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX'.split(',')
# add nygc here when the liftover vcf is ready
sources='pangenie,giggles,pangenie_panel,giggles_panel'.split(',')

wildcard_constraints:
    chr='|'.join(chromosomes),
    sample='|'.join(samples),
    source='|'.join(sources)

### extract sample from nygc phased panel ###
# extract sample from chromosome-wise vcf
rule nygc_extract_sample_chr_vcf:
    input:
        config['path_to_nygc']+'/1kGP_high_coverage_Illumina.{chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz'
    output:
        temp('results/{callset}/callset-comparison/vcfs/nygc-{sample}-{chr}.vcf')
    conda:
        '../envs/basic.yml'
    shell:
        'bcftools view --samples {wildcards.sample} {input} > {output}'

# join chromosome-wise vcfs for a sample
rule nygc_merge_chr_vcfs:
    input:
        expand('results/{{callset}}/callset-comparison/vcfs/nygc-{{sample}}-{chr}.vcf', chr=chromosomes)
    output:
        vcf=temp('results/{callset}/callset-comparison/vcfs/nygc-{sample}.vcf.gz'),
        index=temp('results/{callset}/callset-comparison/vcfs/nygc-{sample}.vcf.gz.tbi')
    conda:
        '../envs/basic.yml'
    shell:
        '''
        bcftools concat -o {output.vcf} -Oz {input}
        tabix -p vcf {output.vcf}
        '''

# extract sample from pangenie genotyped vcf
rule pangenie_create_sample_vcf:
    input:
        config['path_to_pangenie_gts']+'/pangenie-{sample}_genotyping_biallelic.vcf.gz'
    output:
        vcf=temp('results/{callset}/callset-comparison/vcfs/pangenie-{sample}.vcf.gz'),
        index=temp('results/{callset}/callset-comparison/vcfs/pangenie-{sample}.vcf.gz.tbi')
    conda:
        '../envs/basic.yml'
    shell:
        '''
        bcftools view --min-ac 1 {input} | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        '''

# extract sample from giggles genotyped vcf
rule giggles_create_sample_vcf:
    input:
        'results/{callset}/genotypes/{sample}-biallelic.vcf.gz'
    output:
        vcf=temp('results/{callset}/callset-comparison/vcfs/giggles-{sample}.vcf.gz'),
        index=temp('results/{callset}/callset-comparison/vcfs/giggles-{sample}.vcf.gz.tbi')
    conda:
        '../envs/basic.yml'
    shell:
        '''
        bcftools view --min-ac 1 {input} | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        '''

# extract sample vcf from pangenie panel (only for HG01258)
rule pangenie_panel_extract_sample:
    input:
        config['path_to_pangenie_panel_vcf']
    output:
        vcf=temp('results/{callset}/callset-comparison/vcfs/pangenie_panel-HG01258.vcf.gz'),
        index=temp('results/{callset}/callset-comparison/vcfs/pangenie_panel-HG01258.vcf.gz.tbi')
    conda:
        '../envs/basic.yml'
    shell:
        '''
        bcftools view --samples HG01258 {input} | bcftools view --min-ac 1 | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        '''

# extract vcf from giggles panel (only for HG01258)
rule giggles_panel_extract_sample:
    input:
        'results/{callset}/panel/giggles-ready_biallelic.vcf.gz'
    output:
        vcf=temp('results/{callset}/callset-comparison/vcfs/giggles_panel-HG01258.vcf.gz'),
        index=temp('results/{callset}/callset-comparison/vcfs/giggles_panel-HG01258.vcf.gz.tbi')
    conda:
        '../envs/basic.yml'
    shell:
        '''
        bcftools view --samples HG01258 {input} | bcftools view --min-ac 1 | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        '''

# basic counting stats of the vcfs
# set-pass also filters out X and Y chromosomes
rule vcf_stats:
    input:
        'results/{callset}/callset-comparison/vcfs/{source}-{sample}.vcf.gz'
    output:
        "results/{callset}/callset-comparison/vcfs/{source}-{sample}.stats"
    conda:
        "../envs/basic.yml"
    shell:
        "zcat {input} | python scripts/set-pass.py | bcftools view -f PASS --min-af {min_af} --max-af {max_af} | python scripts/vcf_stats.py > {output}"

# add tag of variant types
rule add_tags:
    input:
        vcf='results/{callset}/callset-comparison/vcfs/{source}-{sample}.vcf.gz',
        ref_index=config['path_to_hprc_mg_ref']+'.fai'
    output:
        temp('results/{callset}/callset-comparison/vcfs/{source}-{sample}-{min_af}-{max_af}-tagged.vcf')
    conda:
        "../envs/basic.yml"
    shell:
        "zcat {input.vcf} | python scripts/set-pass.py | bcftools view -f PASS --min-af {wildcards.min_af} --max-af {wildcards.max_af} | python scripts/add-svtags.py | python scripts/truvari-add-contig-len.py {input.ref_index} > {output}"


# intersecting the vcf files and creating upset plots for HG01258 (sample in graph)
rule intersect_HG01258_vcfs:
    input:
        expand('results/{{callset}}/callset-comparison/vcfs/{source}-HG01258-{{min_af}}-{{max_af}}-tagged.vcf', source=sources)
    output:
        tsv="results/{callset}/callset-comparison/in-graph/HG01258-{min_af}-{max_af}/intersection.tsv",
        vcf="results/{callset}/callset-comparison/in-graph/HG01258-{min_af}-{max_af}/intersection.vcf",
        pdf="results/{callset}/callset-comparison/in-graph/HG01258-{min_af}-{max_af}/intersection.pdf"
    conda:
        "../envs/basic.yml"
    log:
        intersect="results/{callset}/callset-comparison/in-graph/HG01258-{min_af}-{max_af}/intersection.log"
    params:
        names=sources
    resources:
        mem_total_mb=4000,
        runtime_hrs=24
    shell:
        '''
        python scripts/intersect-callsets.py intersect -c {input} -n {params.names} -t {output.tsv} -v {output.vcf} -p {output.pdf} &> {log.intersect}
        '''

rule plot_intersect_HG01258:
    input:
        tsv="results/{callset}/callset-comparison/in-graph/HG01258-{min_af}-{max_af}/intersection.tsv"
    output:
        plot="results/{callset}/callset-comparison/in-graph/HG01258-{min_af}-{max_af}/intersection-upsetplot.pdf"
    conda:
        "../envs/basic.yml"
    log:
        plot="results/{callset}/callset-comparison/in-graph/HG01258-{min_af}-{max_af}/plotting.log"
    params:
        columns=["in_" + s for s in sources]
    resources:
        mem_total_mb=4000,
        runtime_hrs=2
    shell:
        '''
        python scripts/plot-comparison-upset.py -t {input.tsv} -o {output.plot} -n {params.columns} &> {log.plot}
        '''

# intersecting the vcf files and creating upset plots for samples not in the graph
rule intersect_outsample_vcfs:
    input:
        expand('results/{{callset}}/callset-comparison/vcfs/{source}-{{sample}}-{{min_af}}-{{max_af}}-tagged.vcf', source=['pangenie', 'giggles'])
    output:
        tsv="results/{callset}/callset-comparison/out-graph/{sample}-{min_af}-{max_af}/intersection.tsv",
        vcf="results/{callset}/callset-comparison/out-graph/{sample}-{min_af}-{max_af}/intersection.vcf",
        pdf="results/{callset}/callset-comparison/out-graph/{sample}-{min_af}-{max_af}/intersection.pdf"
    conda:
        "../envs/basic.yml"
    log:
        intersect="results/{callset}/callset-comparison/out-graph/{sample}-{min_af}-{max_af}/intersection.log"
    params:
        names=['pangenie', 'giggles']
    resources:
        mem_total_mb=4000,
        runtime_hrs=24
    shell:
        '''
        python scripts/intersect-callsets.py intersect -c {input} -n {params.names} -t {output.tsv} -v {output.vcf} -p {output.pdf} &> {log.intersect}
        '''

rule plot_intersect_outsample:
    input:
        tsv="results/{callset}/callset-comparison/out-graph/{sample}-{min_af}-{max_af}/intersection.tsv"
    output:
        plot="results/{callset}/callset-comparison/out-graph/{sample}-{min_af}-{max_af}/intersection-upsetplot.pdf"
    conda:
        "../envs/basic.yml"
    log:
        plot="results/{callset}/callset-comparison/out-graph/{sample}-{min_af}-{max_af}/plotting.log"
    params:
        columns=["in_" + s for s in ['pangenie', 'giggles']]
    resources:
        mem_total_mb=4000,
        runtime_hrs=2
    shell:
        '''
        python scripts/plot-comparison-upset.py -t {input.tsv} -o {output.plot} -n {params.columns} &> {log.plot}
        '''


### truvari comparisons ###
# compare giggles genotypes to pangenie genotypes
rule truvari_outsample_compare:
    input:
        call='results/{callset}/callset-comparison/vcfs/giggles-{sample}-{min_af}-{max_af}-tagged.vcf.gz',
        base='results/{callset}/callset-comparison/vcfs/{source}-{sample}-{min_af}-{max_af}-tagged.vcf.gz'
    output:
        'results/{callset}/truvari-comparison/out-graph/{sample}-{source}-{min_af}-{max_af}/summary.json',
        'results/{callset}/truvari-comparison/out-graph/{sample}-{source}-{min_af}-{max_af}/params.json',
        'results/{callset}/truvari-comparison/out-graph/{sample}-{source}-{min_af}-{max_af}/tp-base.vcf.gz',
        'results/{callset}/truvari-comparison/out-graph/{sample}-{source}-{min_af}-{max_af}/tp-comp.vcf.gz',
        'results/{callset}/truvari-comparison/out-graph/{sample}-{source}-{min_af}-{max_af}/fp.vcf.gz',
        'results/{callset}/truvari-comparison/out-graph/{sample}-{source}-{min_af}-{max_af}/fn.vcf.gz'
    params:
        folder='results/{callset}/truvari-comparison/out-graph/{sample}-{source}-{min_af}-{max_af}',
        tmp='results/{callset}/truvari-comparison/out-graph/{sample}-{source}-{min_af}-{max_af}/tmp'
    conda:
        '../envs/truvari.yml'
    resources:
        mem_total_mb=4000,
        runtime_hrs=2
    shell:
        '''
        truvari bench -b {input.base} -c {input.call} -o {params.tmp}
        mv {params.tmp}/* {params.folder}/
        rm -r {params.tmp}
        '''

# compare giggles genotypes to pangenie genotypes
rule truvari_outsample_compare_lenient:
    input:
        call='results/{callset}/callset-comparison/vcfs/giggles-{sample}-{min_af}-{max_af}-tagged.vcf.gz',
        base='results/{callset}/callset-comparison/vcfs/{source}-{sample}-{min_af}-{max_af}-tagged.vcf.gz',
        ref='/gpfs/project/projects/medbioinf/users/spani/files/fasta/HPRC/reference/chm13v2.0_maskedY_rCRS.fa.gz'
    output:
        'results/{callset}/truvari-comparison-lenient/out-graph/{sample}-{source}-{min_af}-{max_af}/summary.json',
        'results/{callset}/truvari-comparison-lenient/out-graph/{sample}-{source}-{min_af}-{max_af}/params.json',
        'results/{callset}/truvari-comparison-lenient/out-graph/{sample}-{source}-{min_af}-{max_af}/tp-base.vcf.gz',
        'results/{callset}/truvari-comparison-lenient/out-graph/{sample}-{source}-{min_af}-{max_af}/tp-comp.vcf.gz',
        'results/{callset}/truvari-comparison-lenient/out-graph/{sample}-{source}-{min_af}-{max_af}/fp.vcf.gz',
        'results/{callset}/truvari-comparison-lenient/out-graph/{sample}-{source}-{min_af}-{max_af}/fn.vcf.gz'
    params:
        folder='results/{callset}/truvari-comparison-lenient/out-graph/{sample}-{source}-{min_af}-{max_af}',
        tmp='results/{callset}/truvari-comparison-lenient/out-graph/{sample}-{source}-{min_af}-{max_af}/tmp'
    conda:
        '../envs/truvari.yml'
    resources:
        mem_total_mb=4000,
        runtime_hrs=2
    shell:
        '''
        truvari bench -b {input.base} -c {input.call} -o {params.tmp} -f {input.ref} --pick multi -r 2000 --no-ref a -C 2000
        mv {params.tmp}/* {params.folder}/
        rm -r {params.tmp}
        '''


# compare giggles genotypes to pangenie panel (only for HG01258)
rule truvari_HG01258_compare:
    input:
        call='results/{callset}/callset-comparison/vcfs/giggles-HG01258-{min_af}-{max_af}-tagged.vcf.gz',
        base='results/{callset}/callset-comparison/vcfs/{source}-HG01258-{min_af}-{max_af}-tagged.vcf.gz'
    output:
        'results/{callset}/truvari-comparison/in-graph/HG01258-{source}-{min_af}-{max_af}/summary.json',
        'results/{callset}/truvari-comparison/in-graph/HG01258-{source}-{min_af}-{max_af}/params.json',
        'results/{callset}/truvari-comparison/in-graph/HG01258-{source}-{min_af}-{max_af}/tp-base.vcf.gz',
        'results/{callset}/truvari-comparison/in-graph/HG01258-{source}-{min_af}-{max_af}/tp-comp.vcf.gz',
        'results/{callset}/truvari-comparison/in-graph/HG01258-{source}-{min_af}-{max_af}/fp.vcf.gz',
        'results/{callset}/truvari-comparison/in-graph/HG01258-{source}-{min_af}-{max_af}/fn.vcf.gz'
    params:
        folder='results/{callset}/truvari-comparison/in-graph/HG01258-{source}-{min_af}-{max_af}/',
        tmp='results/{callset}/truvari-comparison/in-graph/HG01258-{source}-{min_af}-{max_af}/tmp'
    conda:
        '../envs/truvari.yml'
    resources:
        mem_total_mb=4000,
        runtime_hrs=2
    shell:
        '''
        truvari bench -b {input.base} -c {input.call} -o {params.tmp}
        mv {params.tmp}/* {params.folder}/
        rm -r {params.tmp}
        '''
"""