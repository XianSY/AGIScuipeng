process FUSION_COMPUTE_WEIGHTS{
    tag "$meta.id"
    label 'process_high'

    container "${ 'https://singularity/fusion-gemma' }"
    
    errorStrategy { task.exitStatus in 1..152 ? 'retry' : 'ignore' }


    input:
        tuple val(meta),path(pheno)
        path(snp_binary)

    output:
	path("results"), emit : fusion_weights

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: "--model top1,lasso,enet,blup"

        """
        mkdir -p tmp 
        mkdir -p results
        ln -s . output
        ls ${snp_binary}/*.bed | while read id ; do 
	bfile_name=\$(basename \$id _snp.txt.bed)
	Rscript /opt/fusion/fusion_twas-master/FUSION.compute_weights.R \\
         --bfile ${snp_binary}/\${bfile_name}_snp.txt \\
         --tmp ./tmp/data \\
         --out results/\$bfile_name \\
         --pheno ${pheno}/\${bfile_name}.tx \\
         --noclean TRUE \\
         --PATH_gemma /opt/fusion/gemma_software/gemma-0.98.1-linux-static \\
         $args || { continue; };done >> FUSION.compute_weights.log
       
        [ -f ${pheno}.wgt.RDat ] || touch testwgt.RDat

        """
}
