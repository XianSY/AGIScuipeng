process MENOYPTE2VCF{
    tag "$meta.id"
    label 'process_medium'

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }
    maxRetries 10

    container "${'https://singularity/ewas21'}"

    input:
    tuple val(meta),path(menotype_matrix)

    output:
    tuple val(meta),path("*.SMP.vcf"),emit:SMP_vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cat ${menotype_matrix} | grep -v '##' | awk '
    NR == 1 {
      # 修改第一行的前9列
      printf "##fileformat=VCFv4.2\\n";
      printf "#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT";
      for (i = 5; i <= NF; i++) {
        printf "\\t%s", \$i;
      }
      printf "\\n";
      next;
    }
    {
      # 打印其他行
      printf "%s\\t%s\\t%s:%s\\tC\\tT\\t100\\tPASS\\t.\\tGT", \$1, \$4, \$1, \$4;
      for (i = 5; i <= NF; i++) {
        if (\$i == "MM") printf "\\t0/0";
        else if(\$i == "MU") printf "\\t0/1";
        else if(\$i == "UU") printf "\\t1/1";
        else printf "\\t./.";
      }
      printf "\\n";
    }' > ${menotype_matrix}.SMP.vcf
    
    """
}
