process COMPUTE_ALLETE_RATE{
    tag "$meta.id"
    label "process_low"

    container "${ 'https://singularity/fusion' }"

    input:
        tuple val(meta),path(cis_eqtl)

    output:
        tuple val(meta),path("*.frq"),emit:allele_info

    script:
        """
            plink --vcf ${cis_eqtl} --freq --out allele_info

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                plink: \$(plink --version ) | sed -e "s/PLINK //g" 
            END_VERSION
        
        """

}
