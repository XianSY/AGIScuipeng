process OSCAPCA{
    tag "$meta.id"
    label 'process_medium'

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }
    maxRetries 10

    //装载osca容器
    container "${'https://singularity/osca'}"
    
    input:
    tuple val(meta),path(bod),path(oii),path(opi)
    output:
    tuple val(meta),path("*.eigenval"),path("*.eigenvec"),emit:vec
    when:
    task.ext.when == null || task.ext.when

    script:
    """
    bodname=\$(basename $bod .bod)
    osca --befile \$bodname --pca 10 --out \$bodname
    
    """


}
