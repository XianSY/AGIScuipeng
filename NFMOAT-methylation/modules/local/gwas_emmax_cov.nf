process GWAS_EMMAX_COV {
    tag "${meta.id}"
    label 'process_high'
    // Define input channel

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }

    maxRetries 10 

    container "${'https://singularity/admixture'}"

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'ignore' } 
    maxRetries 10

    input:
    tuple val(meta), path(pfile), path(log), path(tfam)

    // Define output files
    output:
    tuple val(meta), path("*.cov") , emit: cov
    tuple val(meta), path("*.cverror") , emit: cverr     

    when:
    task.ext.when == null || task.ext.when

    script:

    """
        knum=\$(grep "CV error" *K*.log|sed 's/.*=//g;s/): /\\t/g'|sort -gk1|awk '{a=prev-\$2;if(NR>1 && a<0){print k;exit;};prev=\$2;k=\$1;}')
        samplenum=\$(cat $tfam|wc -l)
        if [ "\$samplenum" -le "\$knum" ]; then
            ((knum=samplenum))
        fi
        cat $tfam|awk '{print \$1"\\t"\$1"\\t"1}'|paste - *.\$knum.Q|sed 's/ /\\t/g' > ${meta.id}.emmax.cov
        grep "CV error" *K*.log|sed 's/.*=//g;s/): /\t/g'|sort -gk1|awk '{print \$1,\$2}' > ${meta.id}.emmax.cverror

    """

    
}
