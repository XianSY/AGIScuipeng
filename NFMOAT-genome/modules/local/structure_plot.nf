process STRUCTURE_PLOT{
    tag "${meta.id}"
    label 'process_single'
    // Define input channel

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }

    maxRetries 10 

    container "${'https://singularity/emmax-anno'}"

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'ignore' } 
    maxRetries 10

    input:
    tuple val(meta),path(cverr)

    output:
    tuple val(meta),path("*.png")

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    structure_plot.R \\
    $cverr \\
    >${cverr}.Rplotout
    
    """



}
