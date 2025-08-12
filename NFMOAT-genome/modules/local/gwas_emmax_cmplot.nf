process GWAS_EMMAX_CMPLOT {
    tag "${meta.id}"
    label 'process_single'
    // Define input channel

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }

    maxRetries 10 

    container "${'https://singularity/emmax-anno'}"

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'ignore' } 
    maxRetries 10

    input:
    tuple val(meta), val(pfid), path(qk)

    // Define output files
    output:
    tuple val(meta), path("*.lst") ,emit: signif
    tuple val(meta), path("*.png") ,emit: jpg
    //tuple val(meta), path("*.log") ,emit: log

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    gwas_emmax_cmplot.R \\
    $qk \\
    > ${pfid}.Rout


    """

    
}
