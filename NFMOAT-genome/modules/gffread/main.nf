process GFFREAD{
    input:
    gff

    output:
    path "*.gtf",emit:gtf
    
    script:
    def agrs = task.ext.args ?:''
    def prefix = task.ext.prefix ?:'${gff.baseName}'
    """
    gffread \\
        $gff \\
        $agrs \\
        -o ${prefix}.gtf
        
    """
}