process TWAS_fun{
    tag "${meta.id}"
    label "process_high"

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }

    maxRetries 10
    
    container "${ 'https://singularity/gapit'}"


    input:
        tuple val(meta),path(expression),path(gene_GM),path(phenotype)

    output:
        tuple val(meta),path("*.pdf"),path("*.csv"),path("*")

    script:
    """
        twas.r \
            --expression  ${expression} \
            --gene_GM  ${gene_GM} \
            --phenotype ${phenotype}
    """

}
