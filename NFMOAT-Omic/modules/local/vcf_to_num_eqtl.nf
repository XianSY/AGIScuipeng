process VCF_TO_NUM_EQTL{
    tag "$meta.id"
    label "process_low"

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }

    maxRetries 10

    container "${ 'https://singularity/bedtools2' }"

    input:
        tuple val(meta),path(gt_vcf)

    output:
        tuple val(meta),path("*.num.txt"),path("snps_pos.txt"),emit:vcf_num
    script:
        """
	    plink2 --vcf ${gt_vcf} --export A --out ${gt_vcf}.num;
	    


	    /opt/conda/envs/samtools/bin/awk '{printf "%s", \$2; for(i=7;i<=NF;i++) printf "\\t%s", \$i; print ""}' \
		${gt_vcf}.num.raw > ${gt_vcf}.num.raw.filter

	   /opt/conda/envs/samtools/bin/awk '
    {
        for(i=1;i<=NF;i++){
            a[NR,i]=\$i
        }
        if(NF>max_nf)
            max_nf=NF
    }
    END{
        for(i=1;i<=max_nf;i++){
            for(j=1;j<=NR;j++){
                printf "%s",a[j,i]
                if(j<NR) printf "\\t"
            }
            printf "\\n"
        }
    }' ${gt_vcf}.num.raw.filter > ${gt_vcf}.num.txt
            sed -i 's/#//g' ${gt_vcf}.num.txt

            zcat ${gt_vcf} | awk '!/^##/ {print \$3"_"\$4"\t"\$1"\t"\$2}' > snps_pos.txt
            sed -i 's/#//g' snps_pos.txt 
           
        """
}
