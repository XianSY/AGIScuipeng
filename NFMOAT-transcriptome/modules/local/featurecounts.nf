process FEATURECOUNTS {
    tag "$meta.group"
    label 'process_medium'

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }

    maxRetries 10

    conda "bioconda::subread=2.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/subread:2.0.1--hed695b0_0' :
        'biocontainers/subread:2.0.1--hed695b0_0' }"

    input:
    tuple val(meta),path(annotation)
    path bams  

    output:
    tuple val(meta), path("*featureCounts.txt")        , emit: featCounts
    tuple val(meta), path("*featureCounts.txt.summary"), emit: summary
    tuple val(meta), path("new_counts.txt")            , emit: counts
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def paired_end = meta.single_end ? '' : '-p'

    def strandedness = 0
    if (meta.strandedness == 'forward') {
        strandedness = 1
    } else if (meta.strandedness == 'reverse') {
        strandedness = 2
    }
    """
    featureCounts \\
        $args \\
        $paired_end \\
        -T $task.cpus \\
        -a $annotation \\
        -s $strandedness \\
        -o results.featureCounts.txt \\
        $bams
    sed -i '1 d' results.featureCounts.txt
    sed -i 's/gene://g' results.featureCounts.txt
    sed -i 's/.markdup.sorted.bam//g' results.featureCounts.txt
    sed -i 's/.sorted.bam//g' results.featureCounts.txt
    sed -i 's/.bam//g' results.featureCounts.txt
    awk 'BEGIN{ OFS=","} {printf "%s", \$1; for (i=7; i<=NF; i++) printf " %s", \$i; print ""}' results.featureCounts.txt > new_counts.txt
   sed -i 's/gene://g' new_counts.txt 
   sed -i 's/.markdup.sorted.bam//g' new_counts.txt
   sed -i 's/.bam//g' new_counts.txt 
   cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        subread: \$( echo \$(featureCounts -v 2>&1) | sed -e "s/featureCounts v//g")
    END_VERSIONS
    """
}
