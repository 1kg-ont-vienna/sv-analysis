#compress vcf
rule compress_vcf:
    input:
        "results/{filename}.vcf"
    output:
        vcf="results/{filename}.vcf.gz",
        tbi="results/{filename}.vcf.gz.tbi"
    shell:
        """
        bgzip -c {input} > {output.vcf}
        tabix -p vcf {output.vcf}
        """