process GWAS_EMMAX_ASSOCIATE {
    tag "${meta.id}"
    label 'process_single'
    // Define input channel
   
    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }

    maxRetries 10 

    container "${'https://singularity/emmax'}"

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' } 
    maxRetries 10

    input:
    tuple val(meta), path(tped), path(tfam), path(nosex), path(log), path(kin), path(cov), path(phy), val(whichphy), val(phycol), val(korkq)

    // Define output files
    output:
    tuple val(meta), val("${meta.id}_${whichphy}_${korkq}"), path("*.ps") ,emit: assoc

    when:
    task.ext.when == null || task.ext.when

    script:
        def covarities = cov ? "-c ${cov}" : ''
        def prefix = task.ext.prefix ?: "${meta.id}_${whichphy}_${korkq}"
        def phycolnum=phycol+1
    
    """
    echo ${prefix}
    emmainname=\$(echo $tped| sed 's/\\.tped\$//g')
    awk '{print \$1}' $tfam |xargs -i grep {} $phy| awk -F ',' '{print \$1"\\t"\$1"\\t"\$${phycolnum}}' > ${prefix}.phy

    emmax-intel64 \\
        -t \$emmainname \\
        -o ${prefix}_emmax \\
        -p ${prefix}.phy \\
        -k $kin

    """
    // \\
    //$covarities

    
}
