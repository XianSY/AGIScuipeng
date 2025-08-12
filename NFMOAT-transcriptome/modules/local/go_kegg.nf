process GO_KEGG{
    tag "$meta.id"
    label 'process_medium'

    container "${ 'https://singularity/clusterprofilerwithr' }"

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }
    
    maxRetries 10

    input:
        tuple val(meta),path(gene),val(species)

    output:
        tuple val(meta),path("*csv"),emit:go_kegg

    script:
        """
            go_kegg.r \\
                --input_gene $gene \\
                --species $species
        
        """
}
