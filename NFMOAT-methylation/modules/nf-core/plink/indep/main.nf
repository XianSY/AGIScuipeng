
process PLINK_INDEP {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink:1.90b6.21--h779adbc_1':
        'biocontainers/plink:1.90b6.21--h779adbc_1' }"

    input:
    tuple val(meta), path(vcf), path(tbi),val(chrnum)

    output:
    tuple val(meta), path("*.prune.in")                    , emit: prunein
    tuple val(meta), path("*.prune.out")                   , emit: pruneout
    tuple val(meta), path("*.nosex")                       , emit: prunenosex
    tuple val(meta), path("*.log")                         , emit: prunelog
    path "versions.yml"                                    , emit: versions

    //when:
    //task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    //    --bed ${bed}  \\
    //    --bim ${bim}  \\
    //    --fam ${fam}  \\

    """
    chrn=\$(seq  -s, 1 $chrnum)
    plink \\
        --vcf $vcf \\
        --threads $task.cpus \\
        --indep-pairwise ${task.ext.window_size} ${task.ext.variant_count} ${task.ext.variance_inflation_factor} \\
        --out $prefix \\
        --chr \$chrn \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink: \$(echo \$(plink --version) | sed 's/^PLINK v//;s/64.*//')
    END_VERSIONS
    """


}
