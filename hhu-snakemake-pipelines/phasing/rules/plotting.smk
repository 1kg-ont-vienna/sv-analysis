rule list_trios:
    input:
        expand('results/trio-comparison/{sample}/pairwise.tsv', sample=family_sample_list)
    output:
        'results/trio-comparison/pairwise-list.tsv'
    resources:
        runtime_hrs=0,
        runtime_min=10
    run:
        f = open(output[0], 'w')
        for name in input:
            print(name, file=f)
        f.close()

rule plot_trios:
    input:
        'results/trio-comparison/pairwise-list.tsv'
    output:
        'results/plots/trio-ser.svg',
        'results/plots/trio-ser.pdf'
    params:
        'results/plots/'
    conda:
        '../envs/basic.yaml'
    shell:
        'python scripts/plot-ser-trios.py -tsvs {input} -output {params}'


rule list_non_trios:
    input:
        expand('results/no-trio-comparison/{sample}/pairwise.tsv', sample=samples)
    output:
        'results/no-trio-comparison/pairwise-list.tsv'
    resources:
        runtime_hrs=0,
        runtime_min=10
    run:
        f = open(output[0], 'w')
        for name in input:
            print(name, file=f)
        f.close()

rule plot_non_trios:
    input:
        tsv='results/no-trio-comparison/pairwise-list.tsv',
        meta='resources/igsr_sample_data.tsv'
    output:
        'results/plots/no-trio-ser.svg',
        'results/plots/no-trio-ser.pdf'
    params:
        'results/plots/'
    conda:
        '../envs/basic.yaml'
    shell:
        'python scripts/plot-ser-boxplot.py -tsvs {input.tsv} -meta {input.meta} -output {params}'   

        