process FAMDB{
    tag "$meta.id"
    label "process_low"

    container "${ 'https://singularity/repeatmasker' }"

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }

    maxRetries 10

    input:
        val(meta)

    output:
        tuple val(meta),path("*.embl"),path("*.fa"),emit : ch_famdb

    script:
        def dfam_h5 = params.dfam_h5 ? "-i ${params.dfam_h5}" :"-i /opt/conda/envs/repeatmasker/share/RepeatMasker/Libraries/Dfam.h5"
        """
            /opt/conda/envs/repeatmasker/share/RepeatMasker/famdb.py \\
                $dfam_h5 \\
                families \\
                -f embl  \\
                -a \\
                -d ${params.species} > ${params.species}.embl

            /opt/conda/envs/repeatmasker/share/RepeatMasker/util/buildRMLibFromEMBL.pl ${params.species}.embl > ${params.species}.fa

        """
}
