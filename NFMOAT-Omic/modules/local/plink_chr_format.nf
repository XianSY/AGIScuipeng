process PLINK_CHR_FORMAT{
    tag "$meta.id"
    label "process_low"

    container "${ 'https://singularity/fusion' }"
    
    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }

    maxRetries 10

    input:
        tuple val(meta),path(vcf)
    
    output:
        path 'LDREF',emit:ldref

    script:
        """
            mkdir LDREF

            for chr in {1..${meta.chrome_num}};
            do plink --vcf $vcf --allow-extra-chr --chr \$chr --make-bed --out ./LDREF/popu_snp.\${chr}; done
        
        """
}
