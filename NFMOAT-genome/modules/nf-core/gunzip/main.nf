process GUNZIP {
    input:
    tuple val(meta),path(gzfile)

    output:
    tuple val(meta),path("$gunzip"),emit:gunzip
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?:''
    gunzip = gzfile.toString() - '.gz'

    """
    gzip \\
        -cd \\
        $args \\
        $gzfile \\
        > $gunzip
    
    """
}