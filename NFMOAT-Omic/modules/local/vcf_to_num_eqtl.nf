process VCF_TO_NUM_EQTL{
    tag "$meta.id"
    label "process_low"

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }

    maxRetries 10

    container "${ 'https://singularity/fusion' }"

    input:
        tuple val(meta),path(gt_vcf)

    output:
        tuple val(meta),path("*.num.txt"),path("snps_pos.txt"),emit:vcf_num
    script:
        """
            zcat ${gt_vcf} \
            | awk 'BEGIN{OFS="\\t"} \
                !/^##/ { \
                printf "%s", \$3; \
                for(i=10; i<=NF; i++){ \
                    gt=\$i;\
                    gsub("\\\\.|\\\\./\\\\.","NA",gt); \
                    gsub("^0[\\\\/|]0\$","0",gt); \
                    gsub("^(0[\\\\/|]1|1[\\\\/|]0)\$","1",gt); \
                    gsub("^1[\\\\/|]1\$","2",gt); \
                    printf "%s%s",OFS,gt; \
                    } \
                    printf "\\n"; \
                    }' > ${gt_vcf}.num.txt
            sed -i 's/#//g' ${gt_vcf}.num.txt

            zcat ${gt_vcf} \
                | awk '!/^##/ {print \$3"\t"\$1"\t"\$2}' \
                > snps_pos.txt
            sed -i 's/#//g' snps_pos.txt 
           
        """
}
