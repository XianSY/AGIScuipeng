process GET_ENRICEMENG_GENE{
    tag "enrichment"
    label 'process_single'

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }
    
    maxRetries 10

    conda "${moduleDir}/environment.yml"
    container "${ 'https://singularity/deseq2' }"

    input:
        path(gene_file)

    output:
        path("*txt"),emit:get_gene
    
    script:
        """
            awk -F',' '{print \$1}' ${gene_file} > \$(basename ${gene_file} .csv).txt 
        
        """


}
