process PHENOTYPE_DEAL{
    tag "$meta.id"
    label "process_medium"

    container "${ 'https://singularity/dgeal' }"

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }

    maxRetries 10

    input:
        tuple val(meta),path(expression)

    output:
        path("expression"),emit:expression
	//path 'gene_expression',emit:express_dir
    script:
     	"""
            dgel.r \\
                --expression $expression
            mkdir -p  expression && mv *.tx expression
        """

}
