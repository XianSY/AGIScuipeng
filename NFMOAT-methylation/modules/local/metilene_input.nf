process METILENE_INPUT{
    tag "$meta.id"
    label 'process_medium'

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }
    maxRetries 10

    //装载metilene容器
    container "${'https://singularity/bedtools'}"

    input:
    tuple val(meta),path(group)
    path 'Bed'
    
    output:
    path("*.bedgraph"),emit:bedgraph_matix
    tuple val(meta),path("*.group.bedgraph"),emit:group_bedgraph
    path "versions.yml"                 , emit: version
 
    script:
    """
    bedgraphname=\$(ls -l  Bed1 | awk -F'->' '{print \$2}' | awk -F/ '{print \$NF}' | awk -F'_' '{print \$1}')
    bedtools unionbedg -i \$(ls Bed* | sort) -filler NA -header -names \$(ls -l  Bed* | sort |  awk -F'->' '{print \$2}' | awk -F/ '{print \$NF}' | awk -F'_' '{print \$3}' | xargs) > \$bedgraphname.bedgraph
    bedtools unionbedg -i \$(ls Bed* | sort) -filler NA -header -names \$(ls -l Bed* |sort | awk -F'->' '{print \$2}' | awk -F"/" '{print \$NF}' | awk -F'_' '{print \$3}' | while read var; do awk -v var="\$var" '{a[\$1]=\$2} END {for (i in a) if (i == var) print  a[i]}' $group; done |xargs) > \$bedgraphname.group.bedgraph
     cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bedtools: \$(echo \$(bedtools -v 2>&1))
    END_VERSIONS
   """
}
