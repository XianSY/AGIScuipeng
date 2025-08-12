process BETA_MATRIX_TO_MEQTL{
    tag "$meta.id"
    label "process_low"

    container "${ 'https://singularity/eQTL' }"

    input: 
        tuple val(meta),path(beta_matrix)

    output:
        tuple val(meta),path("smp_beta.txt"),path("smp_pos.txt"),emit:smp_meQTL

    script:
    """
        beta_matrix_to_meQTL.r \\
            --beta_matrix  $beta_matrix    
    """
}
