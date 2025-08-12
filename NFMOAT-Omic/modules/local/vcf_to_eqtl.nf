process VCF_TO_EQTL{
    tag "$meta.id"
    label "process_medium"

    container "${ 'https://singularity/eQTL' }"

    input:
        tuple val(meta),path(vcf_num)

    output:
        tuple val(meta),path("*.txt"),emit:vcf_num
    script:
        """
       	  vcf_transposition.r \\
		--vcf_num $vcf_num

        """
}
