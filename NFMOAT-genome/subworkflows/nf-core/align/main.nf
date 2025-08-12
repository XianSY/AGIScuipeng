include { BWAMEM2_INDEX                 } from '../../../modules/nf-core/bwamem2/index/main'
include { BWAMEM2_MEM                   } from '../../../modules/nf-core/bwamem2/mem/main'  
include { MERGE_FASTQ                   } from '../../../modules/local/merge_fastq'



workflow ALIGN {

    take:
    ch_fastqc_clean
    ch_fasta

    main:
    ch_merge_fastq_files = ch_fastqc_clean
        //ch_merge_fastq_files.view()i
    MERGE_FASTQ (
        ch_merge_fastq_files
    )
    ch_map_fastq = MERGE_FASTQ.out.reads
    //ch_map_fastq.view()
    //def refname = file(fastaPath).baseName
    def bwaidx_meta='fasta'
    ch_fasta_genome_withmeta=Channel.of(bwaidx_meta).combine(ch_fasta).collect().map{it}
    ch_fasta_genome_withoutmeta=ch_fasta
    //ch_fasta_genome_withmeta=tuple(bwaidx_meta, file(params.fasta))
    //ch_fasta_genome_withoutmeta=file(params.fasta)
    //print(ch_fasta_genome_withmeta)
    //print(ch_fasta_genome_withoutmeta)
    
    BWAMEM2_INDEX (
        ch_fasta_genome_withmeta
    )
    //ch_map_fastq.view()
    ch_tomap_fastq=ch_map_fastq.map{
        row -> 
            def sampleid = [:]
            sampleid.id=row[0]
            [sampleid,row[1]]
    }
    ch_bwamen_index = BWAMEM2_INDEX.out.index
    BWAMEM2_MEM (
        ch_tomap_fastq,
        ch_bwamen_index,
        true
    )
    //BWAMEM2_MEM.out.bam.view()
      
    emit:

    bwamem2_mem_bam          = BWAMEM2_MEM.out.bam
}

