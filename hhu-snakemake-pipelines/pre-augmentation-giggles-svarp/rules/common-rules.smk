#compress vcf
rule compress_vcf:
    input:
        "result/svarp-giggles/{filename}.vcf"
    output:
        vcf="result/svarp-giggles/{filename}.vcf.gz",
        tbi="result/svarp-giggles/{filename}.vcf.gz.tbi"
    shell:
        """
        bgzip -c {input} > {output.vcf}
        tabix -p vcf {output.vcf}
        """