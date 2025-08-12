process GGPLOTSTRUCTURE{
    tag "$meta.id"
    label "process_low"

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }
    
    maxRetries 10

    container "${'https://singularity/population'}"

    input:
        tuple val(meta),path(structure),path(group)
    output:
        tuple val(meta),path("*.stat.gz"),path(".png") , emit: lddecay

    script:
        """
            ggplotstructure.r \\
                --input_q_dir  $structure \\
                --input_group_file $group
                    
        """
}