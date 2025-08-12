process SAMTOOLS_FAIDX {
    tag "$fasta"
    label 'process_single'


    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }
    maxRetries 10

    conda "bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0' :
        'biocontainers/samtools:1.17--h00cdaf9_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path ("*.fai"), emit: fai
    tuple val(meta), path ("*.gzi"), emit: gzi, optional: true
    tuple val(meta), path ("*len"), path ("genome_len_kary"), emit: kary
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    samtools \\
        faidx \\
        $args \\
        $fasta

    cut -f 1,2 ${fasta}.fai > chr_len 
    awk -vFS="\t" -vOFS="\t" '{print "chr","-",\$1,"C"NR,"0",\$2,"chr"NR}' chr_len > genome_len_kary     

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch ${fasta}.fai
    cat <<-END_VERSIONS > versions.yml

    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
