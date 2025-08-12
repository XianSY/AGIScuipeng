process GET_CHROMSOME_LST {
    tag "$fasta"
    label 'process_single'

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }

    maxRetries 10 

    input:
    path(fasta)
    //tuple val(meta2), path(fai)

    output:
    path ("*.list") , emit: lst
    //path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
        grep ">" ${fasta}|sed 's/>//g' > allseq.list
        grep ">" ${fasta}|sed 's/>//g'|awk '{if(/^Chr|^[0-9]|^[0-9](2)|^NC_/){print \$1}}' > chrom.list
        grep ">" ${fasta}|sed 's/>//g'|awk '{if(!/^Chr|^[0-9]|^[0-9](2)|^NC_/){print \$1}}' > ctg.list

    """
}
