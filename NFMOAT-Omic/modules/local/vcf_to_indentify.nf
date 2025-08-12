process VCF_TO_INDENTIFY{
    tag "$meta.id"
    label "process_low"

    container "${ 'https://singularity/bedtools' }"

     errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }

    maxRetries 10

    input:
        tuple val(meta),path(vcf)

    output:
        tuple val(meta),path("*.GT.vcf.gz"),emit:gt_vcf
    script:
        """
            bcftools annotate -x FORMAT $vcf -Oz -o ${vcf}.GT.vcf.gz
        
        """
}
