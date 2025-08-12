process MERGE_WEIGHTS_INFO{
    tag "$meta.id"
    label "process_medium"

    container "${ 'https://singularity/fusion' }"

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }

    maxRetries 10 

    input:
        tuple val(meta),path(WEIGHT)
        tuple val(meta),path(genes_bed)

    output:
        tuple val(meta),path("gene_wgt.pos"), emit:gene_wgt
	path (results),emit:WEIGHT

    script:
        """
        
	ls  results/*.wgt.RDat | awk  '{print "'\$(pwd)/'" \$1}'  > gene_weight.txt
	
        awk -F '[/\\.]' '{print \$0"\t"\$(NF-2)}' gene_weight.txt > weight_info.pos
        
         merge_weights_info.R \\
                      --wgt_pos  weight_info.pos \\
                       --bed      $genes_bed
                
        """

}
