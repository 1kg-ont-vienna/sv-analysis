# SV analysis

Long-read sequencing and structural variant characterization in 1,019 samples from the 1000 Genomes Project

## Data collection

The data is hosted at the International Genome Sample Resource (IGSR) in the directory [https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/).

## Genome alignments

The alignment pipeline used to align to GRCh38, CHM13 and a prebuilt human genome graph is detailed below. The [reference subdirectory](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/reference/) at IGSR contains the used linear and graph-based reference genomes.

### Base-calling

Guppy 6.2.1 was used for basecalling.

`guppy_basecaller --disable_pings -x 'cuda:all' -c 'dna_r9.4.1_450bps_sup_prom.cfg' -i fast5 -s fastq --trim_adapters false --do_read_splitting false --compress_fastq false`

### Adapter trimming

Porechop 0.2.4 was used to trim adapters and split reads on internal adapters.

`porechop-runner.py -i in.fasta -o out.fasta --format fasta`

### Linear reference alignments (GRCh38 and CHM13)

Minimap2 v2.26 was used to map the ONT reads.

`minimap2 -ax map-ont --rmq=yes --MD --cs -L -R "@RG\tID:${RG}\tLB:${ID}\tPL:ONT\tPU:${PU}\tSM:${ID}" hg38.fa in.fasta`

Samtools v1.16.1 was used to sort alignments and convert to cram.

`samtools sort -m3G --reference hg38.fa -O cram -o ${RG}.hg38.cram in.sam`

Multiple ONT runs for the same sample were tagged using different read-groups and merged using samtools v1.16.1.

`samtools merge --reference hg38.fa -O cram -o ${ID}.hg38.cram ${RGs}.hg38.cram`

### Graph genome alignments

The pangenome graph built with minigraph for HPRC year-1 samples (https://doi.org/10.5281/zenodo.6983934) was used to map the ONT reads with minigraph.

`minigraph --vc -cx lr chm13-90c.r518.gfa.gz ${ID}.fasta | bgzip > ${ID}.gaf.gz`

## Data Reuse policy

This README relates to all data associated with the collaborative effort analyzing the long-read sequencing data of the 1019 members of the 1KG-ONT panel. These data were generated at the Institute of Molecular Pathology (Vienna, Austria) with funds provided by Boehringer-Ingelheim. Extracted DNA was obtained from the Coriell Institute for Medical Research and was consented for full public release of genomic data. Please see Coriell (https://www.coriell.org) for more information on specific cell lines.

All data in or under the directory at [https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/) is part of this Project.

If using this data please acknowledge that: These data were generated at the Institute of Molecular Pathology (Vienna, Austria) with funds provided by Boehringer-Ingelheim.

Continuing the philosophy of the 1000 Genomes Project, the data producer of this Project will release the Project data quickly, prior to publication, in the expectation that they will be valuable for many researchers. In keeping with Fort Lauderdale principles, data users may use the data for many studies, but are expected to allow the data producers to make the first presentations and to publish the first paper with global analyses of the data.


### Global analyses of Project data

The Project plans to publish global analyses of the sequence data and quality, structural variants, STRs, microsatellites. haplotypes, LD patterns, population genetic phenomena such as population comparisons, mutation rates, signals of selection, and functional annotations, as well as analyses of regions of general interest. Talks, posters, and papers on all such analyses are to be published first by the Project team. When the first major Project paper on these analyses is published, then researchers inside and outside the Project are free to present and publish using the Project data for these and other analyses.

### Samples for methods development and evaluation

A subset of 1KG_ONT_VIENNA genomes overlap with samples analyzed by the Human Genome Structural Variation Consortium, namely HG00268, HG00513, HG00731, HG02554, HG02953, NA12878, NA19129, NA19238, NA19331, NA19347. These samples can be used for methods development and evaluation prior to a publication by the project team.

Researchers who have questions about whether they may make presentations or submit papers using Project data can contact us.

For any queries related to the data or the data reuse policy, please contact us at 1kg.ont@imp.ac.at
