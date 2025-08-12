process ANNO_TABLE {
    tag "${meta.id}"
    label 'process_single'

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }
    
    maxRetries 10

    //conda "${moduleDir}/environment.yml"
    container "${'https://singularity/annovar-with-perl'}"

    input:
    tuple val(meta), path(vcf),path(tbi), path(fa), path(fa2), path(gff), path(list),path(anno)

    output:
    tuple val(meta), path("*.exonic_variant_function"), emit: annoevf
    tuple val(meta), path("*.variant_function"), emit: annovf
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    
    """
    annotate_variation.pl \\
        -geneanno \\
        -dbtype refGene \\
        -out $prefix \\
        -build $meta.id\\
        $anno \\
        .

    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        table_annovar: \$(echo \$(table_annovar.pl version 2>&1) )
    END_VERSIONS
    """

    
    //BUILD_DOC=`find -L ./ -name "*.list" | sed 's/\\${meta.id}_refGene.list\$//'`
    /*
    """
    table_annovar.pl \\
        $vcf \\
        . \\
        -buildver $meta.id \\
        -out ${prefix}\\
        -protocol refGene \\
        -operation g \\
        -nastring . \\
        -vcfinput
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        table_annovar: \$(echo \$(table_annovar.pl version 2>&1) )
    END_VERSIONS
    """
    */
}
