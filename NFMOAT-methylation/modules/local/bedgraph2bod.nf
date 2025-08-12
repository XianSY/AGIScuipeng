process BEDGRAPH2BOD{
    tag "$meta.id"
    label 'process_high'

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }
    maxRetries 10

    //装载metilene容器
    container "${'https://singularity/metilene-2.8'}"

    input:
    //tuple val(meta),path(chg),path(cpg),path(chh)
    tuple val(meta),path(bedgraph)    
  

    output:
    tuple val(meta),path("*.txt"),emit:bod
    //tuple val(meta),path("CPG*.txt"),emit:cpgbod
    //tuple val(meta),path("CHH*.txt"),emit:chhbod

    when:
    task.ext.when == null || task.ext.when


    script:
    """
    awk '{a=\$1;b=\$2;\$2="";\$1=a"_"b;\$3="";print \$0}' $bedgraph | awk '{
  for (i = 1; i <= NF; i++) {
    a[i, NR] = \$i
  }
  max_col = (NF > max_col ? NF : max_col)
  max_row = NR
} 
END {
  for (i = 1; i <= max_col; i++) {
    for (j = 1; j <= max_row; j++) {
      printf "%s%s", a[i, j], (j == max_row ? "\\n" : " ")
    }
  }
}' | awk '{print \$1,\$0}' > ${bedgraph}.txt
   
    
    """

}
