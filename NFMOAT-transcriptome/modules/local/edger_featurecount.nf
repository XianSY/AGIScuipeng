process EDGER_COUNTS{
    tag "$meta.group"
    label 'process_medium'

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }

    maxRetries 10

    conda "${moduleDir}/environment.yml"
    container "${ 'https://singularity/edger' }"

    input:
    tuple val(meta),path(counts),path(compare),path(group)

    output:
    tuple val(meta),path("*.csv"),path("*.png"),path("*.pdf"),emit:deseq

    when:
    task.ext.when == null || task.ext.when
    
    script:
    """
    edgeR.r \\
        --counts_file $counts \\
        --group_file $group \\
        --compare_file $compare \\
        --outdir ./

    """


}
