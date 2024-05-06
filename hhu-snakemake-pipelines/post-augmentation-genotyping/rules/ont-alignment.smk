# extract fasta from cram
rule cram_to_fasta:
    input:
        cram=config['path_to_cram']+'/{sample}.hg38.cram',
        ref=config['reference_directory']+'/1KG_ONT_VIENNA_hg38.fa'
    output:
        fasta=temp("result/data/fasta/{sample}.fasta.gz"),
        fai=temp("result/data/fasta/{sample}.fasta.gz.fai"),
        gzi=temp("result/data/fasta/{sample}.fasta.gz.gzi")
    params:
        ref_dir=config['reference_directory']
    conda:
        "../envs/basic.yml"
    resources:
        runtime_hrs=20,
        mem_total_mb=20000
    priority: 1
    shell:
        '''
        seq_cache_populate.pl -root {params.ref_dir} {input.ref}
        export REF_PATH={params.ref_dir}/%2s/%2s/%s:http://www.ebi.ac.uk/ena/cram/md5/%s
        export REF_CACHE={params.ref_dir}/%2s/%2s/%s
        samtools fasta {input.cram} | bgzip -c > {output.fasta}
        samtools faidx {output.fasta}
        '''

# run minigraph alignment
rule minigraph_alignment:
    input:
        fasta='results/data/fasta/{sample}.fasta.gz',
        fai='results/data/fasta/{sample}.fasta.gz.fai',
        gzi='results/data/fasta/{sample}.fasta.gz.gzi',
        ref='results/{callset}/bubble_calling_and_tagging/{callset}.tagged.gfa',
        path_to_minigraph=config['path_to_minigraph']
    output:
        temp('results/{callset}/ont-alignments/{sample}.gaf')
    log:
        'results/{callset}/ont-alignments/{sample}.alignment.log'
    resources:
        runtime_hrs=36,
        runtime_min=0,
        mem_total_mb=96000
    threads: 8
    shell:
        '{input.path_to_minigraph} --vc -cx lr {input.ref} {input.fasta} -t {threads} > {output}'


# gaf sorting
rule gaftools_sort:
    input:
        gaf='results/{callset}/ont-alignments/{sample}.gaf',
        gfa='results/{callset}/bubble_calling_and_tagging/{callset}.tagged.gfa'
    output:
        sorted_gaf='results/{callset}/ont-alignments/{sample}.sorted.gaf.gz',
        index='results/{callset}/ont-alignments/{sample}.sorted.gaf.gz.gai'
    log:
        'results/{callset}/ont-alignments/{sample}.sorting.log'
    resources:
        runtime_hrs=54,
        runtime_min=0,
        mem_total_mb=20000
    shell:
        '''
        set +u
        source ~/.bashrc
        conda activate gaftools-env
        set -u
        gaftools sort --bgzip --outgaf {output.sorted_gaf} {input.gaf} {input.gfa} 2> {log}
        '''