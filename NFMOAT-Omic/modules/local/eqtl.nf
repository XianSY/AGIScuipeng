process EQTL_PROCESS{
    tag "$meta.id"
    label "process_high"

    container "${ 'https://singularity/matrixeqtl' }"
     errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }

    maxRetries 10

    input:
        tuple val(meta),path(vcf_num),path(vcf_pos)
        tuple val(meta),path(expression),path(gene_info)

    output:
        tuple val(meta),path("*")

    when:
        task.ext.when == null || take.ext.when

    script:
        """ 
	    ouputname=\$(basename $vcf_num .txt)
            eQTL.r \\
                --gene_info $gene_info \\
                --snp_pos $vcf_pos \\
                --snp_vcf_num $vcf_num \\
                --gene_expression $expression \\
                --output_name \$ouputname
        	
        """


}
