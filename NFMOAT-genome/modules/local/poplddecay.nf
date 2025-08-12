process POPLDDECAY{
    tag "$meta.id"
    label "process_high"

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }
    
    maxRetries 10

    container "${'https://singularity/population'}"

    input:
        tuple val(meta),path(group),path(vcf)
    output:
        tuple val(meta),path("*.stat.gz"),path("*.pdf") , emit: lddecay

    script:
        """
            awk -F ',' 'NR>1 {print \$1 >> \$2".txt"}' $group

            ls *.txt | while read id ; do \\
                /opt/PopLDdecay/bin/PopLDdecay -InVCF ${vcf} -OutStat \${id}.stat.gz -SubPop \${id} ; done
            
            ls *.stat.gz | awk -F '.' '{print \$0"\\t"\$1}' >>population.list

            perl /opt/PopLDdecay/bin/Plot_MultiPop.pl -inList population.list --output LDdecay -bin1 1 -bin2 100
                    
        """
}
