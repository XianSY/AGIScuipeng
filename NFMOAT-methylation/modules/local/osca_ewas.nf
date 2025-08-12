process OSCAEWAS{
    tag "$meta.id"
    label 'process_medium'

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }
    maxRetries 10

    //装载osca容器
    container "${'https://singularity/osca'}"

    input:
    tuple path(phenotype),val(meta),path(bod),path(oii),path(opi),path(pca)

    output:
    tuple val(meta),path("*.linear"),emit:ewas_results
    tuple val(meta),path("*.log"),emit:ewas_log   
 
    when:

    script:
    """
    typename=\$(basename $bod .bod)
    osca \\
        --befile \$typename \\
        --pheno $phenotype \\
        --qcovar $pca \\
        --linear \\
        --out \$typename.${phenotype}
    """
}
