process SEPARATE_SNP_FOR_GENE{
    tag "$meta.id"
    label "process_low"

    container "${ 'https://singularity/fusion' }"

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }

    maxRetries 10

    input:
          tuple val(meta),path(eQTL)

    output:
          path "snp",emit:snp_dir
    when:
          task.ext.when == null || task.ext.when

    script:
        """
            mkdir snp
            cat ${eQTL} | \
            awk 'NR!=1 {print \$2}' | \
            sort | uniq | \
            while read id; do awk -v id="\$id" '\$2==id {print \$1}' ${eQTL} > \${id}_snp.txt ; done
            mv *_snp.txt snp 
        """
          
}
