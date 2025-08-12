process FILTER_PASS_VARIANT {
    tag "$meta.id"
    label 'process_single'
    
    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }
    
    maxRetries 10

    container "${'https://singularity/filter-pass-variant'}"
//    conda "${moduleDir}/environment.yml"
//    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
//        'https://depot.galaxyproject.org/singularity/gatk4:4.4.0.0--py36hdfd78af_0':
//        'biocontainers/gatk4:4.4.0.0--py36hdfd78af_0' }"


    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"


    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Fillter pass variant] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    gzip -cd ${vcf} | perl -ne 'print if /^#/ or /PASS/' | bgzip -c > ${prefix}.vcf.gz
    
    """
}
