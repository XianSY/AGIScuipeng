process VCF_TO_NUMBER{
    tag "$meta.id"
    label "process_low"

    container "${ 'https://singularity/eQTL' }"

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }

    maxRetries 10

    input:
        tuple val(meta),path(vcf)

    output:
        tuple val(meta),path("*.raw"),path("*.txt"),path("*.pos"),emit:vcf_num

    script:
        """
            plink --vcf $vcf --recode A --allow-extra-chr   --out \$(basename ${vcf} .gz).num

            less -S \$(basename ${vcf} .gz).num.raw | gawk '{printf \$2 " "; for (i=7; i<=NF; i++) printf \$i " "; print ""}' > ${vcf}_num.txt

	    bgzip -dc $vcf | grep -v "^##" | awk '{ print \$3, \$1, \$2 }' > ${vcf}.pos
            sed -i 's/#//g' ${vcf}.pos
            sed -i 's/_T//g' ${vcf}_num.txt
            sed -i 's/_G//g' ${vcf}_num.txt
            sed -i 's/_C//g' ${vcf}_num.txt
            sed -i 's/_A//g' ${vcf}_num.txt
        """
}
