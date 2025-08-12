process GETCHRLENGTH{
    tag "$reference_index"
    label "process_low"

    container "${'https://singularity/histogram'}"

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }
    maxRetries 10


    input:
        path(reference_index)

    output:
        path("contigs.txt"),emit:contig

    script:
        """
            cut -f1,2 ${reference_index} | awk '{print "##contig=<ID="\$1",length="\$2">";}' > contigs.txt
            echo "##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">" >> contigs.txt
        """

}
