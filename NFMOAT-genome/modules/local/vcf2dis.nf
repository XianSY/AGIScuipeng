process VCF2DIS{
    tag "$meta.id"
    label "process_low"

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }
    
    maxRetries 10

    container "${'https://singularity/vcf2dis'}"

    input:
        tuple val(meta),path(vcf)

    output:
        tuple val(meta),path("*.nwk"),emit:nwk

    script:
        """
            VCF2Dis -InPut ${vcf} -OutPut ${vcf}.mat
        
        """


}
