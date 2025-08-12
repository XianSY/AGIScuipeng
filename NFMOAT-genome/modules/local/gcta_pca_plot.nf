process GCTA_PCA_PLOT{
    tag "${meta.id}"
    label 'process_single'
    // Define input channel

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }
    
    maxRetries 10

    container "${'https://singularity/emmax-anno'}"

    input:
    tuple val(meta),path(eigenvec)

    output:
    tuple val(meta),path("*.png")

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    gcta_pca_plot.R \\
    $eigenvec \\
    >${eigenvec}.Rplotout
    
    """



}
