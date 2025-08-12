process MAKE_LD{
    tag "$meta.id"
    label "process_medium"

    container "${ 'https://singularity/dgeal' }"

    when:
        task.ext.when == null | task.ext.when

    input:
        tuple val(meta),path(allele_info),path(gwas_results)

    output:
        tuple val(meta),path("*.sumstats"),emit:sumstats

    script:
        """
            make_LD.R \\
                --gwas_resluts $gwas_results \\
                --allele_info $allele_info 


            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
            R version : \$(Rscript --version) 
            END_VERSION
        """
}
