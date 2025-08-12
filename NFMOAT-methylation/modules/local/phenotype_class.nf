process PHENOTYPE_CLASS{
     tag "$meta.id"
    label 'process_medium'

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }
    maxRetries 10

    container "${'https://singularity/histogram'}"

    input:
    tuple val(meta),path(phenotype)

    output:
    path("*.txt"),emit:pheno

    script:
    """
    phenotype_class.r \\
       --input_phenotype $phenotype

    """
}
