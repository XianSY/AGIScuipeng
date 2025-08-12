process PLOTKINSHIP{
    tag "$meta.id"
    label "process_low"

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }
    
    maxRetries 10

    container "${'https://singularity/population'}"

    input:
        tuple val(meta),path(kif)
    output:
        tuple val(meta),path("*.stat.gz"),path(".png") , emit: lddecay

    script:
        """
            plotkinship.r \\
                --input_Kif ${kif}
                    
        """
}