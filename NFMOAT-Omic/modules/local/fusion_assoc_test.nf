process FUSION_ASSOC_TEST{
    tag "$meta.id"
    label "process_medium"

    container "${ 'https://singularity/fusion-gemma' }"

    
    input:
        path(ldchr)
 	tuple val(metadata),path(sumstats)
        path(weight_dir)
        tuple val(meta),path(gene_wgt)
    output:
        path("twas_results*")    
    when:
        task.ext.when == null || tash.ext.when

    script:
        """
               for chr in {1..${metadata.chrome_num}};
               do  
                 Rscript /opt/fusion/fusion_twas-master/FUSION.assoc_test.R  \\
                --sumstats $sumstats \\
                --weights $gene_wgt \\
                --weights_dir  / \\
                --ref_ld_chr ${ldchr}/popu_snp. \\
                --chr \$chr \\
                --out twas_results_\${chr}.dat
	      done
        """
}
