process VCF_TO_BIALLELIC{
    tag "$meta.id"
    label "process_low"

    container "${ 'https://singularity/bedtools' }"

    input:
        tuple val(meta),path(vcf)

    output:
        tuple val(meta),path("*.biallelic.vcf.gz"),emit:biallelic_vcf
    script:
        """
            bcftools \
                view -m2 -M2 -v snps $vcf \
                -Oz -o ${vcf}.biallelic.vcf.gz
        
        """
}