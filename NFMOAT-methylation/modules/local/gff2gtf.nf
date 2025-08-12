process GFF2GTF{
    tag "${gff}"
    label "process_low"

    container "${ 'https://singularity/bedtools' }"

    input:
        path(gff)
    output:
        path("*.gtf"),emit: gtf
    script:
    """
        gffread -T $gff > ${gff}.gtf
    """
}
