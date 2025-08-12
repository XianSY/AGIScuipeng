process GFF2GTF{
    tag "$meta.id"
    label "process_low"

    container "${ 'https://singularity/bedtools' }"

    input:
        tuple val(meta),path(gff)
    output:
        tuple val(meta),path("*.gtf"),emit: gtf
    script:
    """
        gffread -T $gff > ${gff}.gtf
    """
}
