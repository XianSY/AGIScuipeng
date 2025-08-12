
process PLINK_MAKEBED {
    tag "$meta.id"
    label 'process_low'

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }

    maxRetries 10 

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink:1.90b6.21--h779adbc_1':
        'biocontainers/plink:1.90b6.21--h779adbc_1' }"

    input:
    tuple val(meta), path(ped), path(map)
    //val(window_size)
    //val(variant_count)
    //val(variance_inflation_factor)

    output:
    tuple val(meta), path("*.bed"), path("*.bim"), path("*.fam") ,path("*.nosex") ,path("*.log"), emit: bed
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    pedname=\$(echo $ped| sed 's/\\.ped\$//g')
    plink \\
        --file \$pedname\\
        --threads $task.cpus \\
        --make-bed --allow-extra-chr --noweb \\
        --out \$pedname \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink: \$(echo \$(plink --version) | sed 's/^PLINK v//;s/64.*//')
    END_VERSIONS
    """
}
