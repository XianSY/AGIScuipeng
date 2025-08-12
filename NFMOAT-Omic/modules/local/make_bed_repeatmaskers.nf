process MAKE_BED_REPEATMASKER{
    tag "$meta.id"
    label "process_low"

    container "${'https://singularity/bedtools'}"

     errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }

    maxRetries 10

    input:
        tuple val(meta),path(vcf)
    output:
        tuple val(meta),path("*.bed"),emit:bed
    script:

          """
          cat $vcf | grep '/' | awk -v OFS="\t"  '{print \$5,\$6,\$7,\$11}' >\$(basename $vcf .out).bed
	  sed -i '1d' \$(basename $vcf .out).bed       
 
          """
}
