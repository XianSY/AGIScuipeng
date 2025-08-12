process FILTER_BED {
    tag "$meta.id"
    label 'process_high'

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }
    maxRetries 10

    container "${'https://singularity/histogram'}"

    input:
    tuple val(meta),path(bedgraph)

    output:
    tuple val(meta),path("*.filter.bed") , emit:filter_bed

    script:

    """
    filter_bedgraph.r \\
        --bedgraph_file $bedgraph
    
    """


}
