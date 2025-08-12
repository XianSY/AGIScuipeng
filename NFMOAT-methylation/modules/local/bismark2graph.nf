process BISMARK2BEDGRAPH {
    tag "$meta.id" 
    label 'process_high'

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }
    maxRetries 10

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bismark:0.24.0--hdfd78af_0' :
        'biocontainers/bismark:0.24.0--hdfd78af_0' }"

    input:
    tuple val(meta),path(methylation_calls),path(cx_report)
    
    output:
    tuple val(meta),path("*.bedgraph.gz"),   emit: bedgraph
    tuple val(meta),path("*.bismark.cov.gz"),emit: bedgraph_cov
    path "versions.yml"                     ,emit: version
     
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    echo $methylation_calls
    output_name=\$(basename $methylation_calls .txt).bedgraph
    echo output_name
    bismark2bedGraph \\
         --CX $cx_report \\
         -o \$output_name \\
         --buffer_size 70G \\
         ${args}           \\
         $methylation_calls

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bismark: \$(echo \$(bismark -v 2>&1) | sed 's/^.*Bismark Version: v//; s/Copyright.*\$//')
    END_VERSIONS
    """



}
