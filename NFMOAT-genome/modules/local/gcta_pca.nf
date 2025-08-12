process GCTA_PCA {
    tag "${meta.id}"
    label 'process_single'
    // Define input channel
    
    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }
    
    maxRetries 10

    container "${'https://singularity/gcta'}"

    input:
    tuple val(meta), path(grmnbin), path(grmbin) ,path(grmid) ,path(log)
    // Define output files
    output:
    tuple val(meta), path("*.eigenval") ,path("*.eigenvec") , emit: vec

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    samplenum=\$(cat $grmid|wc -l)
    grmname=\$(echo $grmid|sed 's/\\.grm\\.id//g')
    gcta64 \\
    --grm \$grmname \\
    --pca \$samplenum \\
    --out \$grmname

    """

    
}
