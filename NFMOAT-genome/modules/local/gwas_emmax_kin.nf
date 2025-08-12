process GWAS_EMMAX_KINSHIP {
    tag "${meta.id}"
    label 'process_single'
    // Define input channel

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }

    maxRetries 10 

    container "${'https://singularity/emmax'}"

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'ignore' } 
    maxRetries 10

    input:
    tuple val(meta), path(tped), path(tfam), path(nosex), path(log)

    // Define output files
    output:
    tuple val(meta), path("*.kinf") ,emit: kin

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    emmainname=\$(echo $tped| sed 's/\\.tped\$//g')

    emmax-kin-intel64 \\
        \$emmainname \\
        -v \\
        -d 10 \\
        -o ${meta.id}.BN.kinf


    """

    
}
