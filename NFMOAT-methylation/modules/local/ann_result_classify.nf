process ANN_RESULTS_CLASSIFY{
    tag "$meta.id"
    label 'process_medium'

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }
    maxRetries 10

    container "${'https://singularity/bedtools'}"

    input:
    tuple val(meta), path(ann_results),path(chr_len),path(kary)

    output:
   // path "*.map.out"

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cat $ann_results | awk -F'[:\t|]' '\$5=="upstream_gene_variant"{print \$1,\$2-1,\$2,\$5}' >upstream_gene_variant.bed
    cat $ann_results | awk -F'[:\t|]' '\$5=="intron_variant"{print \$1,\$2-1,\$2,\$5}' >intron_variant.bed
    cat $ann_results | awk -F'[:\t|]' '\$5=="missense_variant"{print \$1,\$2-1,\$2,\$5}' >missense_variant.bed
    cat $ann_results | awk -F'[:\t|]' '\$5=="intergenic_region"{print \$1,\$2-1,\$2,\$5}' >intergenic_region.bed
    cat $ann_results | awk -F'[:\t|]' '\$5=="5_prime_UTR_variant"{print \$1,\$2-1,\$2,\$5}' >5_prime_UTR_variant.bed
    
    grep -v '^[SAPM]' $chr_len > chr_len.txt    

    bedtools makewindows -g chr_len.txt -w 3000000 > len.100k

    ls *.bed | while read id;do bedtools coverage -a len.100k -b \$id | cut -f1,2,3,4 > \$id.map.out; done
    """
}
