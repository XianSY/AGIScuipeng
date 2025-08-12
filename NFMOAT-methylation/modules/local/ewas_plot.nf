process EWAS_PLOT{
    tag "$meta.id"
    label 'process_medium'

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }
    maxRetries 10

    container "${'htpps://singularity/ewas-plot'}"

    input:
    tuple val(meta),path(ewas_result)

    output:
    tuple val(meta),path("*jpg")

    script:
    """
    ewas_plot.r \\
       --input_ewas_results $ewas_result
    
    """

}
