process VCF_TO_INDENTIFY{
    tag "$meta.id"
    label "process_low"

    container "${ 'https://singularity/bedtools2' }"

     errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }

    maxRetries 10

    input:
        tuple val(meta),path(vcf)

    output:
        tuple val(meta),path("*.GT.vcf.gz"),emit:gt_vcf
	tuple val(meta),path("*eigenvec.T"),  emit:eigenvec
    script:
        """
            bcftools annotate -x FORMAT $vcf -Oz -o ${vcf}.GT.vcf.gz
	    plink2 --vcf $vcf --pca 10 --out ${vcf}.pca
	    
	    sed -i 's/#IID/IID/g' ${vcf}.pca.eigenvec
	    
	    awk '{ for (i=1;i<=NF;i++) {
                    a[i,NR]=\$i
                }
                max_nf=NF
                max_nr=NR
                }
            END {
                for (i=1;i<=max_nf;i++) {
                    for (j=1;j<=max_nr;j++) {
                        printf a[i,j]
                        if (j<max_nr) printf "\\t"
                    }
                    printf "\\n"
                }
            }' ${vcf}.pca.eigenvec > ${vcf}.pca.eigenvec.T
        
        """
}
