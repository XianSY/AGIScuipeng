process STURCTURE_VCF_ADDID {
    tag "${meta.id}"
    label 'process_single'
    // Define input channel
  
    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }

    maxRetries 10 
   
    container "${'https://singularity/bedtools'}"

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'ignore' } 
    maxRetries 10

    input:
    tuple val(meta), path(vcf)

    // Define output files
    output:
    tuple val(meta), path("*sort.vcf.gz") ,path("*.tbi"),env(chrnum),emit: vcf
    tuple val(meta), env(allctgnum), emit: ctgnum
    when:
    task.ext.when == null || task.ext.when

    script:

    """
        mkdir -p ./work/tmp
        
        chrnum=\$( awk '{print \$1}' $vcf | sort | uniq | grep -v '#' | wc -l)
	allctgnum=\$( grep -E "^#" $vcf | wc -l)
        awk -v FS='\\t' -v OFS='\\t' '{if(/^#/){print \$0}else{a=\$1;gsub(/[Ctg|Chr]/,"",a);\$3=a":"\$2;print \$0}}' $vcf | \
        bgzip -c > ${meta.id}.vcf.gz

        bcftools sort ${meta.id}.vcf.gz -Oz -o ${meta.id}.sort.vcf.gz
        tabix ${meta.id}.sort.vcf.gz

        echo "chrnum_output=\$chrnum" > ${meta.id}.chrnum.txt
    """

    
}
