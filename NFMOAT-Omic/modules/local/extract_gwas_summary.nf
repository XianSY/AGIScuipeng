process EXTRACT_GWAS_SUMMARY{
    tag "$meta.id"
    label "process_low"

    container "${ 'https://singularity/bedtools' }"

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }

    maxRetries 10


    input:
        tuple val(meta),path(gwas_results),path(vcf)

    output:
        tuple val(meta),path("gwas_results_summary"),emit:gwas_summary

    script:
        """
           cat ${gwas_results} | grep -v "nan" | \
           awk '{print \$1"\t"\$3/\$2}' > \$(basename ${gwas_results} .ps).summary_info

           cat ${gwas_results} | grep -v "nan" | \
           awk '{print \$1}' | sort  > \$(basename ${gwas_results} .ps).snpId

           vcftools --gzvcf ${vcf} --snps \$(basename ${gwas_results} .ps).snpId --recode --stdout \
                | bcftools query -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\n" > \$(basename ${gwas_results} .ps).base_info

           sort -k1,1 \$(basename ${gwas_results} .ps).summary_info > \$(basename ${gwas_results} .ps).sort.summary_info

           sort -k3,3 \$(basename ${gwas_results} .ps).base_info > \$(basename ${gwas_results} .ps).sort.base_info

           join -t \$'\t' -1 3 -2 1 \
               -o 1.3 1.4 1.5 2.2 \
               \$(basename ${gwas_results} .ps).sort.base_info \
               \$(basename ${gwas_results} .ps).sort.summary_info > gwas_results_summary


            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                bedtools: \$(bedtools --version ) | sed -e "s/PLINK //g" 
            END_VERSION
        
        """

}
