process GATK4_ANALYZECOVARIATES{
    tag "$meta.id"
    label 'process_medium'

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }
    maxRetries 10 


    conda "bioconda::gatk4=4.4.0.0"
    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'https://depot.galaxyproject.org/singularity/gatk4:4.4.0.0--py36hdfd78af_0':
    //    'biocontainers/gatk4:4.4.0.0--py36hdfd78af_0' }"
    container "${'https://singularity/gatk4withr'}"

    input:
    tuple val(meta), path(input_before_bqsr), path(input_after_bqsr)

    output:
    tuple val(meta), path("*.pdf")       , emit: pdf
    //path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 3072

    if (!task.memory) {
        log.info '[GATK AnalyzeCovariates] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M" AnalyzeCovariates \\
        -before $input_before_bqsr \\
        -after $input_after_bqsr \\
        -plots ${prefix}.pdf \\
        --tmp-dir . \\
        $args
    """
    /*
    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
    */
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
   
    /*
    """
    touch ${prefix}.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
    */
}
