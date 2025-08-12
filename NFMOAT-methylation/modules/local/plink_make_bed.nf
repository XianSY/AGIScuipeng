process PLINK_MAKEBED{
    tag "$meta.id"
    label 'process_low'

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }

    maxRetries 10 

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink:1.90b6.21--h779adbc_1':
        'biocontainers/plink:1.90b6.21--h779adbc_1' }"

    input:
        tuple val(meta),path(bed), path(bim), path(fam)
    
    output:
         tuple val(meta), path("*.tped"), path("*.tfam"), path("*.nosex"), path("*.log"), emit: emmaxfile
         path "versions.yml"                                    , emit: versions
    
    script:
        """
            name=\$(echo $bed | sed 's/\\.bed\$//g')
            plink --bfile \$name \\
                --recode 12 transpose \\
                --allow-extra-chr \\
                --out \${name}.emmax_in 

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                plink: \$(echo \$(plink --version) | sed 's/^PLINK v//;s/64.*//')
            END_VERSIONS
        
        """
}
