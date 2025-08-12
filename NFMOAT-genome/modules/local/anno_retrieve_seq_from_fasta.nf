process ANNO_RETRIEVE_SEQ {
    tag "${meta.id}"
    label 'process_single'

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }
    
    maxRetries 10

    //conda "${moduleDir}/environment.yml"
    container "${'https://singularity/annovar-with-perl'}"

    input:
    tuple val(meta), path(genepred), path(fasta)

    output:
    tuple val(meta), path("*_refGeneMrna.fa"), emit: fa
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    
    """
    retrieve_seq_from_fasta.pl \\
        --format refGene \\
        --seqfile $fasta \\
        $genepred \\
        --out ${prefix}.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        retrieve_seq_from_fasta: \$(echo \$(retrieve_seq_from_fasta.pl version 2>&1) )
    END_VERSIONS
    """

}
