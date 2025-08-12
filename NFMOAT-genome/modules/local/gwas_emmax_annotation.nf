process GWAS_EMMAX_ANNOTATION {
    tag "${meta.id}"
    label 'process_single'
    // Define input channel

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }

    maxRetries 10 

    container "${'https://singularity/admixture'}"

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'ignore' } 
    maxRetries 10

    input:
    tuple val(meta), path(signif), path(annov)

    // Define output files
    output:
    tuple val(meta), path("*.anno") ,emit: anno

    when:
    task.ext.when == null || task.ext.when

    script:

    """
        less -S $annov| \
        awk -v FS='\\t' -v OFS='\\t' '{if(\$4==a && \$5==b){\$2=c"##"\$2;\$3=d"##"\$3;\$7=e"##"\$7;\$8=f"##"\$8;}else{print g;}a=\$4;b=\$5;c=\$2;d=\$3;e=\$7;f=\$8;g=\$0}END{print g}'|
        tail -n +2|
        awk -F '\\t' '{a=\$3;gsub(/:.*\$/,"",a);b=\$4;gsub(/[Chr|Ctg]/,"",b);gsub(/gene-/,"",a);print b":"\$5"\\t"a"\\t"\$0}'|
        sort -k 1b,1 > ${signif}.annotations_groupbysnp.lst

        sort -k 1b,1 $signif|join -i - ${signif}.annotations_groupbysnp.lst|sort -gk3,1 -gk2,2 -s > ${signif}.anno

    """

    
}
