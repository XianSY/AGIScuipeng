process VCF_ADD_HEADER{
    tag "${meta.id}"
    label "process_low"

    container "${'https://singularity/bedtools'}"

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }
    maxRetries 10


    input:
        tuple val(meta),path(menotype_vcf)
        path(contig)

    output:
        tuple val(meta),path("*.header.vcf"),emit:SMP_vcf

    script:
        """
            output_name=\$(basename ${menotype_vcf} .vcf.gz).header.vcf
            bcftools annotate \\
                -h ${contig} \\
                -Oz -o \$output_name \\
                ${menotype_vcf}
        """

}
