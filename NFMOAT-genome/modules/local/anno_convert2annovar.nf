process ANNO_CONVERT2ANNOVAR {
    tag "${meta.id}"
    label 'process_single'

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }
    
    maxRetries 10

    //conda "${moduleDir}/environment.yml"
    container "${'https://singularity/annovar-with-perl'}"

    input:
    tuple val(meta), path(vcf),path(tbi)

    output:
    tuple val(meta), path("*.var"), emit: anno
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    
    """
    convert2annovar.pl \
        -format vcf4 \
        -allsample \
        -withfreq \
        $vcf > ./${prefix}.var

    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        table_annovar: \$(echo \$(table_annovar.pl version 2>&1) )
    END_VERSIONS
    """
}
