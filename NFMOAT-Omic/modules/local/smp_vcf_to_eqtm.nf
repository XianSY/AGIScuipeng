process VCF_TO_EQTM{
    tag "$meta.id"
    label "process_medium"

    container "${ 'https://singularity/eQTL' }"

    input:
        tuple val(meta),path(vcf_num)

    output:
        tuple val(meta),path("popu_smp.chromesome_t.txt"),path("smps_pos.txt"),emit:vcf_num
    script:
        """
	    smp_vcf_transposition.r \\
		--vcf_num $vcf_num
 
        """
}
