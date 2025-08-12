process EMMAX_VCF_GROUPPING {
    tag "$meta.id"
    label 'process_single'

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }

    maxRetries 10 

    container "${'https://singularity/vcftools-groupping'}"

    input:
    tuple val(meta), val(gsample), path(vcf), path(tbi)

    // Define output files
    output:
    tuple val(meta), path("*.vcf.gz") ,path("*.tbi") ,emit: vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def indv_command = gsample.collect(itpth -> "--indv ${itpth}").join(' ')
    //println("good")
    """
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
    vcftools \\
        --gzvcf $vcf \\
        --recode \\
        --stdout | bgzip -c > ${meta.id}_group.vcf.gz
    
    tabix ${meta.id}_group.vcf.gz

    """
}
