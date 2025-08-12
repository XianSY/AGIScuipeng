process PLOT_HISTORGAM{
    tag "$meta.id"
    label "process_medium"

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }
    maxRetries 10

    container "${'https://singularity/histogram'}"

    input:
    tuple val(meta),path(bedgraph)
    val(cx)
    
    output:
    tuple val(meta),path("*.png")

    script:
    """

    plot_histogram.r \\
        --bedgraph $bedgraph \\
        --cx $cx
    
    """
}
