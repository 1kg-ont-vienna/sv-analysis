configfile: 'config.yaml' 

# tag rGFA chromosome-wise
rule call_rGFA_bubbles:
    input:
        ref=lambda wildcards: config['callsets'][wildcards.callset]['gfa']
    output:
        expand('results/{{callset}}/bubble_calling_and_tagging/{{callset}}-{chr}.gfa', chr=chromosomes),
        expand('results/{{callset}}/bubble_calling_and_tagging/{{callset}}-{chr}.csv', chr=chromosomes)
    params:
        out_dir='results/{callset}/bubble_calling_and_tagging'
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

# concat chromsome-wise tagged GFA
rule concat_tagged_GFA:
    input:
       lambda wildcards: expand('results/{{callset}}/bubble_calling_and_tagging/%s-{chr}.gfa'%(config['callsets'][wildcards.callset]['gfa'].split("/")[-1][:-4]), chr=chromosomes)
    output:
        'results/{callset}/bubble_calling_and_tagging/{callset}.tagged.gfa'
    shell:
        'cat {input} > {output}'