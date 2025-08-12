process CORRECT_INPUT{
    label "process_low"

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }

    container "${'https://singularity/bedtools'}"


    input:
    path(fasta)
    path(gff)

    output:
    path "*cor.fasta" , emit : out_fasta
    path "*cor.gff" , emit : out_gff

    script:
        """
            sed  's/^>Chr/>/' ${fasta} > \$(basename ${fasta} .fasta).cor.fasta
            awk 'BEGIN{FS=OFS="\t"} {sub(/^Chr/, "", \$1); print}' ${gff} > \$(basename ${gff} .gff).cor.gff
        
        """


}
