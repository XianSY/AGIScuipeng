process PLINK_FORMAT{
    tag "$meta.id"
    label "process_high"

    container "${ 'https://singularity/fusion' }"

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }

    maxRetries 10

    input:
        tuple val(meta),path(vcf)
        path(snp_dir) 
    output:
        path "snp_binary",emit:plink_format       
 
    script:
        """
            mkdir -p snp_binary
            ls snp | while read id ; do plink --vcf ${vcf} --extract snp/\$id --allow-extra-chr  --make-bed --out snp_binary/\$id; done


            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                plink:\$(plink --version ) | sed -e "s/PLINK //g")
	    END_VERSIONS
        """    
}
