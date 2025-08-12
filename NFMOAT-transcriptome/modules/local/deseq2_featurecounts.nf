process DESEQ2_COUNTS{
    tag "$meta.group"
    label 'process_medium'

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }
    
    maxRetries 10

    conda "${moduleDir}/environment.yml"
    container "${ 'https://singularity/deseq2' }"

    input:
    tuple val(meta),path(counts),path(compare),path(group)

    output:
    path("*.csv"),emit:deseq_results
    tuple val(meta),path("*.png"),path("*.pdf"),emit:deseq_image

    when:
    task.ext.when == null || task.ext.when
    
    script:
    """
    deseq2_counts.r \\
        --counts_file $counts \\
        --group_file $group \\
        --compare_file $compare \\
        --outdir ./

    """


}
