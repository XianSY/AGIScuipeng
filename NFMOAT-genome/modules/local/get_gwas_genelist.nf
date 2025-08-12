process GET_ENRICEMENG_GENE{
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ 'https://singularity/bedtools' }"

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }
    
    maxRetries 10

    input:
        path(gene_file)

    output:
        path("*txt")
    
    script:
        """
           cat popu_snp_agr09_KQ_emmax.ps_significant.lst.anno | cut -d ' ' -f8 | cut -d ':' -f2 > \$(basename ${gene_file} .lst.anno).txt  
        
        """


}
