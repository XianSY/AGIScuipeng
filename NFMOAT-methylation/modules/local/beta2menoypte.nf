process BETA2MENOYPTE{
    tag "$meta.id"
    label 'process_medium'

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }
    maxRetries 10

    container "${'https://singularity/ewas21'}"
    
    input:
      tuple val(meta),path(beta_matrix)

      output:
       tuple val(meta),path("*menotype"), emit: smp_menoypte
      
      script:

       """
        java -jar /opt/ewas21/ewas2.1.jar \\
                -SMP.convert  \\
                -input  $beta_matrix \\
                -output ${beta_matrix}.menotype
       """
}
