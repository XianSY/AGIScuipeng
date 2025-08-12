
process ADMIXTURE {
    tag "$meta.id"
    label 'process_single'

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'ignore' }
    maxRetries 10

    conda "${moduleDir}/environment.yml"
    container "${'https://singularity/admixture'}"

    input:
    tuple val(meta), path (bed), path(bim), path(fam), val(K)

    output:
    tuple val(meta), path("*.Q")    , emit: ancestry_fractions
    tuple val(meta), path("*.P")    , emit: allele_frequencies
    tuple val(meta), path("*.log")  , emit: cvlog
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"


    """
    admixture \\
        --cv $bed \\
        $K \\
        -j$task.cpus \\
        >> ${meta.id}_K${K}.log
  

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        admixture: \$(echo \$(admixture 2>&1) | head -n 1 | grep -o "ADMIXTURE Version [0-9.]*" | sed 's/ADMIXTURE Version //' )
    END_VERSIONS

    """
}
