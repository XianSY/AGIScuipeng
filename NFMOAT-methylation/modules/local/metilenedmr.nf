process METILENEDMR{
    tag "$meta.id"
    label "process_single"

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }
    maxRetries 10

    //装载metilene容器
    container "${'https://singularity/metilene-2.8'}"

    input:
    tuple val(A),val(B)
    tuple val(meta),path(chg),path(cpg),path(chh)

    output:
    tuple val(meta),path("*.output"),emit:metilene_dmr
    path "versions.yml"                 , emit: version

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    chg_name=(\$(basename $chg .group.bedgraph).metilene.group.bedgraph)
    cp $chg ./\$chg_name 
    cpg_name=(\$(basename $cpg .group.bedgraph).metilene.group.bedgraph)
    cp $cpg ./\$cpg_name
    chh_name=(\$(basename $chh .group.bedgraph).metilene.group.bedgraph)
    cp $chh ./\$chh_name
    metilene \\
         -a $A \\
         -b $B  \\
         -t 8 \\
         \$chg_name > metilene_chg_dmr_$A_vs_$B.output

    metilene \\
         -a $A \\
         -b $B  \\
         -t 8 \\
         \$chh_name > metilene_chh_dmr_$A_vs_$B.output

    metilene \\
         -a $A \\
         -b $B  \\
         -t 8 \\
         \$cpg_name > metilene_cpg_dmr_$A_vs_$B.output

    cat <<-END_VERSIONS > versions.yml
     "${task.process}":
         metilene: \$(echo "0.2-8")
     END_VERSIONS

    """

}
