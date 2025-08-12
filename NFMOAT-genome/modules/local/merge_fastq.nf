process MERGE_FASTQ {
    tag "$sampleID"
    label 'process_single'
    // Define input channel

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }

    maxRetries 10 

 
    input:
    tuple val(sampleID), val(sampleType),val(single_end), path(fastq1), path(fastq2)

    // Define output files
    output:
    
    tuple val(sampleID), path('*.fastq.gz') , emit: reads

    when:
    task.ext.when == null || task.ext.when

    script:
    
    if (sampleType.size() > 1)
    {
        // Merge _1.fq.gz and _2.fq.gz separately when multiple types exist
    """
        cat ${fastq1} > ${sampleID}_1.fastq.gz
        cat ${fastq2} > ${sampleID}_2.fastq.gz
    """
    } else {
        // Create symbolic links for single types
    """
        ln -s ${fastq1} ${sampleID}_1.fastq.gz
        ln -s ${fastq2} ${sampleID}_2.fastq.gz
    """
    }
    
}
