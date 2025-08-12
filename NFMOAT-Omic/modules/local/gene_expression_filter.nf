process GENE_EXPRESSION_FILTER{
    tag "$meta.id"
    label "process_low"

    container "${ 'https://singularity/fusion' }"

    input:
        tuple val(meta),path(gene_info),path(expression)
    
    output:
        tuple val(meta),path("*.tab"),path("gene_info.txt"),emit:gene_expression_eQTL

    script:
        """
            gene_expression_filter.R \\
                --expression  $expression \\
                --gene_info   $gene_info   
        """

}
