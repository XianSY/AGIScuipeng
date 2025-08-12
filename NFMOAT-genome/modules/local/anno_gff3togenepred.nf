process ANNO_GFF3TOGENEPRED {
    tag "${meta.id}"
    label 'process_single'

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }
    
    maxRetries 10


    //conda "${moduleDir}/environment.yml"
    container "${'https://singularity/annovar-with-perl'}"

    input:
    tuple val(meta), path(gff3)

    output:
    tuple val(meta), path("*_refGene.txt"), emit: genep
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    gff3ToGenePred \\
        $gff3 \\
        ${prefix} \\


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gff3ToGenePred: \$(echo \$(gff3ToGenePred version 2>&1) )
    END_VERSIONS
    """

}
