process MAKE_BED{
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
          gunzip -c   $vcf |grep -v "^#" | awk '{
            if(length(\$4) == length(\$5)) {
            print \$1 "\t" \$2 "\t" \$2+1 "\t"\$1":"\$2
            } else {
                print \$1 "\t" \$2 "\t" (\$2 + length(\$4))"\t"\$1":"\$2
            }
        }' > \$(basename $vcf .vcf.gz).bed
          """
}
