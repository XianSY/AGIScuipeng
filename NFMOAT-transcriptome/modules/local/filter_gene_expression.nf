process FILTER_GENE_EXPRESSION{
    tag "$meta.group"
    label 'process_medium'

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }
    
    maxRetries 10

    conda "${moduleDir}/environment.yml"
    container "${ 'https://singularity/deseq2' }"

    input:
    tuple val(meta),path(counts)

    output:
    tuple val(meta),path("filter_gene_expression.txt"),emit:filter_gene_expression

    when:
    task.ext.when == null || task.ext.when
    
    script:
    """
        filter_gene_expression.r \\
            --gene_counts $counts 

    """
}
