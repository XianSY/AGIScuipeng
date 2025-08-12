process BISMARK_METHYLATIONEXTRACTOR {
    tag "$meta.id"
    label 'process_high'

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }
     maxRetries 10


    conda "bioconda::bismark=0.24.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bismark:0.24.0--hdfd78af_0' :
        'biocontainers/bismark:0.24.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam)
    path index

    output:
    tuple val(meta), path("*.bedGraph.gz")         , emit: bedgraph
    tuple val(meta), path("CpG_*.deduplicated.txt"),path("*.CX_report.txt")          , emit: methylation_cpg
    tuple val(meta), path("CHG_*.deduplicated.txt"),path("*.CX_report.txt")          , emit: methylation_chg
    tuple val(meta), path("CHH_*.deduplicated.txt"),path("*.CX_report.txt")          , emit: methylation_chh
    tuple val(meta), path("*.cov.gz")              , emit: coverage
    tuple val(meta), path("*_splitting_report.txt"), emit: report
    tuple val(meta), path("*.M-bias.txt")          , emit: mbias
    tuple val(meta), path("*.cytosine_context_summary.txt") , emit: cytosine_context_summary
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Assign sensible numbers for multicore and buffer_size based on bismark docs
    if(!args.contains('--multicore') && task.cpus >= 6){
        args += " --multicore ${(task.cpus / 2) as int}"
    }
    // Only set buffer_size when there are more than 6.GB of memory available
    if(!args.contains('--buffer_size') && task.memory?.giga > 6){
        args += " --buffer_size ${task.memory.giga - 2}G"
    }

    def seqtype  = meta.single_end ? '-s' : '-p'
    """
    bismark_methylation_extractor \\
        $bam \\
        --comprehensive \\
        --bedGraph \\
        --counts \\
        --report \\
        --CX_context \\
        --cytosine_report \\
        --genome_folder $index \\
        $seqtype \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bismark: \$(echo \$(bismark -v 2>&1) | sed 's/^.*Bismark Version: v//; s/Copyright.*\$//')
    END_VERSIONS
    """
}
