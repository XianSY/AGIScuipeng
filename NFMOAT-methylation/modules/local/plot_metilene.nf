process PLOT_METILENE{
    tag "$meta.id"
    label 'process_medium'

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }
    maxRetries 10

    //装载metilene容器
    container "${'https://singularity/methylseqanno'}"


    input:
    tuple val(meta),path(group_bedgraph),path(gtf)

    output:
    
    tuple val(meta),path("*.csv"),path("*.png"),emit:anno_results


    script:
    """
    awk '{print \$1"\\t"\$2"\\t"\$3}' $group_bedgraph > anno.txt
    ChIPSeeker_anno.r \\
       anno.txt \\
       $gtf
    """

}
