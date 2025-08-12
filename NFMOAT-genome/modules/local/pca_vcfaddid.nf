process PCA_VCF_ADDID {
    tag "${meta.id}"
    label 'process_single'
    // Define input channel

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }

    maxRetries 10 

    container "${'https://singularity/bedtools'}"

    input:
    tuple val(meta), path(vcf), path(tbi)

    // Define output files
    output:
    tuple val(meta), path("*.addID.vcf.gz") ,path("*vcf.gz.tbi"),emit: vcf

    when:
    task.ext.when == null || task.ext.when

    script:

    """
        chrnum=\$(bgzip -dc $vcf | grep -v "^#" | awk '{print \$1}' | grep -E '^[0-9]' | sort | uniq | paste -sd, )
        bcftools view -r \$chrnum $vcf -o filter_${vcf}

        gunzip -c filter_${vcf} | awk 'BEGIN {OFS="\t"} /^#/ {print; next} { \$3 = \$1 ":" \$2; print }' | bgzip -c >\$(basename filter_${vcf} .vcf.gz).addID.vcf.gz
        tabix -p vcf \$(basename filter_${vcf} .vcf.gz).addID.vcf.gz

    """

    
}
