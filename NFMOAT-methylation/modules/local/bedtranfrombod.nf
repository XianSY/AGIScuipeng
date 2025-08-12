process BEDTRANFORMBOD{
    tag "$meta.id"
    label 'process_medium'

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }
    maxRetries 10

    //装载metilene容器
    container "${'https://singularity/osca'}"

    input:
    tuple val(meta),path(chg),path(cpg),path(chh)
    
    output:
    tuple val(meta),path("CHG.bod"),path("CHG.oii"),path("CHG.opi"),emit:CHG_BOD
    tuple val(meta),path("CPG.bod"),path("CPG.oii"),path("CPG.opi"),emit:CPG_BOD
    tuple val(meta),path("CHH.bod"),path("CHH.oii"),path("CHH.opi"),emit:CHH_BOD

    script:
    """
    osca --efile $chg --make-bod --out CHG --methylation-beta
    osca --efile $cpg --make-bod --out CPG --methylation-beta
    osca --efile $chh --make-bod --out CHH --methylation-beta
    """
}
