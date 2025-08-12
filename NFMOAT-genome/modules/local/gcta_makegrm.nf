process GCTA_MAKEGRM {
    tag "${meta.id}"
    label 'process_single'
    // Define input channel

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }
    
    maxRetries 10

    container "${'https://singularity/gcta'}"

    input:
    tuple val(meta), path(ped), path(map), path(bed), path(bim), path(fam) ,path(nosex) ,path(log)

    // Define output files
    output:
    tuple val(meta), path("*grm.N.bin"), path("*.grm.bin") ,path("*.grm.id") ,path("*.log")  ,emit: grm

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    pedname=\$(echo $ped| sed 's/\\.ped\$//g')
    gcta64 \\
    --bfile \$pedname \\
    --make-grm \\
    --out grm_\$pedname


    """

    
}
