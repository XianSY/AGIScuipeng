
process GWAS_EMMAX_PLINK_FMAT {
    tag "$meta.id"
    label 'process_low'

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' } 
    maxRetries 10

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink:1.90b6.21--h779adbc_1':
        'biocontainers/plink:1.90b6.21--h779adbc_1' }"

    input:
    tuple val(meta), path(vcf), path(tbi), val(chrnum), val(allctgnum)
    //val(window_size)
    //val(variant_count)
    //val(variance_inflation_factor)

    output:
    tuple val(meta), path("*.tped"), path("*.tfam"), path("*.nosex"), path("*.log"), emit: emmaxfile
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    chrn=\$(seq -s, 1 $chrnum)
    plink \\
        --vcf $vcf \\
        --recode 12 transpose \\
        --output-missing-genotype 0 \\
        --out emmmax_${prefix} \\
        --allow-extra-chr \\
        --chr \$chrn \\
        $args
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink: \$(echo \$(plink --version) | sed 's/^PLINK v//;s/64.*//')
    END_VERSIONS
    """
}
