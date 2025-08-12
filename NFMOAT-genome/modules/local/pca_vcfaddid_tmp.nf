process PCA_VCF_ADDID {
    tag "${meta.id}"
    label 'process_single'
    // Define input channel
    container "${'https://singularity/admixture'}"

    input:
    tuple val(meta), path(vcf), path(tbi)

    // Define output files
    output:
    tuple val(meta), path("*.vcf.gz") ,path("*.tbi"), env(chrnum) ,emit: vcf

    when:
    task.ext.when == null || task.ext.when

    script:

    """
        bgzip -dc $vcf|\
        awk -v FS='\t' -v OFS='\t' '{if(/^##contig=<ID=/){a=\$0;gsub(/[Ctg|Chr]/,"",a);print a;}else if(/^#/){print \$0;}else{a=\$1;gsub(/[a-zA-Z]+/,"",a);\$1=a;\$3=a":"\$2;print \$0}}'|\
        bgzip -c > ${meta.id}.vcf.gz

        tabix ${meta.id}.vcf.gz


    """

    
}
