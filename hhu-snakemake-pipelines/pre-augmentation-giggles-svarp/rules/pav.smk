###### PAV rules from Arda ######
rule all:
    input:
        'pav_hg38/run.complete',
        'pav_t2t/run.complete'

rule prepare_pav_assembly_table:
    input:
        svtigs_h1=config['sample']+'_svtigs_H1.fa',
        svtigs_h2=config['sample']+'_svtigs_H2.fa',
    output:
        svtigs_h1_hg38=temp('pav_hg38/svtigs_h1.asm.fa'),
        svtigs_h2_hg38=temp('pav_hg38/svtigs_h2.asm.fa'),
        assm_tsv_hg38='pav_hg38/assemblies.tsv',
        svtigs_h1_t2t=temp('pav_t2t/svtigs_h1.asm.fa'),
        svtigs_h2_t2t=temp('pav_t2t/svtigs_h2.asm.fa'),
        assm_tsv_t2t='pav_t2t/assemblies.tsv'
    run:
        import pathlib as pl
        import shutil
        # I don't trust Snakemake keeping a sort-order intact
        input_assemblies = [input.svtigs_h1, input.svtigs_h2]
        output_assemblies_hg38 = [output.svtigs_h1_hg38, output.svtigs_h2_hg38]
        output_assemblies_t2t = [output.svtigs_h1_t2t, output.svtigs_h2_t2t]
        # relative path to working directory for PAV assemblies.tsv
        pav_input_assemblies_hg38 = [fp.replace('pav_hg38/', '', 1) for fp in output_assemblies_hg38]
        pav_input_assemblies_t2t = [fp.replace('pav_t2t/', '', 1) for fp in output_assemblies_t2t]
        names = [pl.Path(fp).name.split('_')[0] for fp in output_assemblies_hg38]
        with open(output.assm_tsv_hg38, 'w') as table:
            _ = table.write('\t'.join(['NAME', 'HAP1', 'HAP2']) + '\n')
            _ = table.write('\t'.join([names[0], pav_input_assemblies_hg38[0], pav_input_assemblies_hg38[1]])  + '\n')
        with open(output.assm_tsv_t2t, 'w') as table:
            _ = table.write('\t'.join(['NAME', 'HAP1', 'HAP2']) + '\n')
            _ = table.write('\t'.join([names[0], pav_input_assemblies_t2t[0], pav_input_assemblies_t2t[1]])  + '\n')
        for infile, outfile in zip(input_assemblies, output_assemblies_hg38):
            shutil.copy(infile, outfile)
        for infile, outfile in zip(input_assemblies, output_assemblies_t2t):
            shutil.copy(infile, outfile)

rule prepare_pav_config:
    input:
        assm_tsv='pav_{ref}/assemblies.tsv',
        ref=config['path_to_ref']+'/1KG_ONT_VIENNA_{ref}.fa',
        ref_idx=config['path_to_ref']+'/1KG_ONT_VIENNA_{ref}.fa.fai'
    output:
        cfg='pav_{ref}/config.json',
        ref=temp('pav_{ref}/data/ref/{ref}.fa'),
        ref_idx=temp('pav_{ref}/data/ref/{ref}.fa.fai'),
    params:
        ref='pav_{ref}/'
    run:
        import json
        import shutil
        pav_cfg = {
            'assembly_table': input.assm_tsv.split('/')[-1],
            'reference': output.ref.replace(params.ref, '', 1)
        }
        with open(output.cfg, 'w') as dump:
            _ = json.dump(pav_cfg, dump, ensure_ascii=True)
        shutil.copy(input.ref, output.ref)
        shutil.copy(input.ref_idx, output.ref_idx)

rule run_pav:
    input:
        'pav_{ref}/config.json',
        'pav_{ref}/svtigs_h1.asm.fa',
        'pav_{ref}/svtigs_h2.asm.fa',
        'pav_{ref}/data/ref/{ref}.fa',
        'pav_{ref}/data/ref/{ref}.fa.fai'
    output:
        chk='pav_{ref}/run.complete'
    params:
        wdir='pav_{ref}/'
    threads: 4
    log:
        pav='pav_{ref}/pav.log'
    benchmark:
        'pav_{ref}/pav_run.rsrc'
    container:
        config['path_to_pav']
    shell:
        'snakemake --verbose --jobs {threads} -d {params.wdir} -s /opt/pav/Snakefile --rerun-incomplete --keep-incomplete --restart-times 0 &> {log.pav} && touch {output.chk}'