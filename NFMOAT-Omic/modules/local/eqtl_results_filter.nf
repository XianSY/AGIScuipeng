process EQTL_RESULTS_FILTER{
    tag "$meta.id"
    label "process_high"

    container "${ 'https://singularity/eQTL' }"
    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }
    
    maxRetries 10

    input:
        tuple val(meta),path(vcf_num),path(vcf_pos)
        tuple val(meta),path(eqtl_results)
        tuple val(meta),path(expression),path(gene_info) 
    output:
        tuple val(meta),path("*")
    
    when:
        task.ext.when == null || take.ext.when
    
    script:
        """
           eQTL_results.r \\
                --input_MatrixeQTL_results $eqtl_results \\
                --snp_info $vcf_pos \\
                --gene_info $gene_info \\
                --output ./
        
        
        """

}
