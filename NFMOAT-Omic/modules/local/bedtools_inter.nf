process BEDTOOLS_INTER{
    tag "$meta.id"
    label "process_low"

    container "${'https://singularity/bedtools'}"

    input:
        tuple val(meta),path(pure_methy),path(inter_file)
    output:
        tuple val(meta),path("*.bed"),emit:pure_methy_file
    script:
        """
            bedtools intersect -a $pure_methy -b $inter_file -v > pure_methy_\$(basename $inter_file .vcf.gz).bed 

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
            bedtools: \$(echo \$(bedtools --version 2>&1))
            END_VERSIONS
        """
}
