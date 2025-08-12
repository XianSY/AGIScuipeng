process GATK4_GATHERVCFS {
    tag "$meta.id"
    label 'process_high'

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }
    maxRetries 10

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.4.0.0--py36hdfd78af_0':
        'biocontainers/gatk4:4.4.0.0--py36hdfd78af_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val(meta), path("*.{tbi,idx}"), emit: tbi
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_list = vcf.collect{"-I ${it}"}.join(' ')

    def avail_mem = 16000
    if (!task.memory) {
        log.info '[GATK GatherVcfs] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        GatherVcfs \\
        $input_list \\
        -O ${prefix}.vcf.gz \\
        $args
    
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        IndexFeatureFile \\
        --input ${prefix}.vcf.gz \\
        --tmp-dir . \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
