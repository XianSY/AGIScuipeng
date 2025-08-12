process MAKE_GENE_INFO{
    tag "$meta.id"
    label "process_low"

    container "${ 'https://singularity/bedtools' }"
    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }

    maxRetries 10 
 
    input:
        tuple val(meta),path(gtf)
    output:
        tuple val(meta),path("genes.bed"),emit:gene_info
    script:
    """
        export HOME=./
        export MPLCONFIGDIR=./
        export PYGTFTK_DIR=./
        mkdir -p \$MPLCONFIGDIR
        convrt_gtf_name=\$(basename ${gtf} .gtf).covrt.gtf
        gtftk convert_ensembl -i ${gtf} | sed -E 's/^(chr|Chr)//' > \$convrt_gtf_name
    
	awk '\$3=="gene" {
        gsub(/[\";]/, "", \$0)
        split(\$10, a, "[;:]")
        count=0
        for(i in a){
          count++
        }
        print a[count] "\\t" \$1 "\\t" \$4-1 "\\t" \$5
    }' \$convrt_gtf_name > genes.bed   
    """
}
