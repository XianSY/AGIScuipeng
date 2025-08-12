process TPM_PCA{
    tag "$meta.group"
    label 'process_medium'

     errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }
    
    maxRetries 10

    conda "${moduleDir}/environment.yml"
    container "${ 'https://singularity/deseq1' }"

    input:
    tuple val(meta),path(counts),path(group)

    output:
    tuple val(meta),path("TPM.csv"),path("gene_info.txt"),path("*.png"),emit:deseq

    when:
    task.ext.when == null || task.ext.when
    
    script:
    """
    pca_plot.r \\
        --counts_file $counts \\
        --group_file $group \\
        --outdir ./


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-deseq2: \$(Rscript -e "library(DESeq2); cat(as.character(packageVersion('DESeq2')))")
    END_VERSIONS
    """


}
