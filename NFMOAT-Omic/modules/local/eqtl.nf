process EQTL_PROCESS{
    tag "$meta.id"
    label "process_high"

    container "${ 'https://singularity/matrixeqtl' }"
     errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }

    maxRetries 10

    input:
        tuple val(meta),path(vcf_num),path(vcf_pos)
        tuple val(meta),path(expression),path(gene_info)
	tuple val(meta),path(cvrt_file)
    output:
        tuple val(meta),path("*.eQTL") , emit:init_eqtl

    when:
        task.ext.when == null || take.ext.when

    script:
        """ 
	    ouputname=\$(basename $vcf_num .vcf.gz.num.txt).eQTL
            eQTL.r \\
                --gene_info $gene_info \\
                --snp_pos $vcf_pos \\
                --snp_vcf_num $vcf_num \\
                --gene_expression $expression \\
                --output_name \$ouputname \\
        	--cvrt_file $cvrt_file
        """


}
