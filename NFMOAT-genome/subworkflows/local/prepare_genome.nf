 include{ GUNZIP as GUNZIP_FASTA         }from '../../modules/nf-core/gunzip'
 include{ GUNZIP as  GUNZIP_GTF          }from '../../modules/nf-core/gunzip'
 include{ GUNZIP as GUNZIP_ADAPTER_FASTA }from '../../modules/nf-core/gunzip'
 include{ GUNZIP }from '../../modules/nf-core/gunzip/main'
 include{ CORRECT_INPUT                  }from '../../modules/local/correct_input'

workflow PREPARE_GENOME {
    take:
    fasta                     //file : genome file .fa
    gtf                       //file : gtf file .gtf
    gff                       //file : gff file .gtf
    adapter_fasta             //file : adapter file .fa

    main:
    
    
    // ungzip the fasta.gz file
    if (fasta.endsWith(".gz")) {
        // ungip fasta.gz file 
        ch_fasta = GUNZIP_FASTA([[:],fasta]).gunzip.map{it[1]}
        
    }else {
        ch_fasta = Channel.value(file(fasta))
        
    }
    
    // ungzip gtf.gz || gff.gz file 
    if ( gtf || gff ) {
        if (gtf) {
            if ( gtf.endsWith(".gz")){
                ch_gtf = GUNZIP_GTF([[:],gtf]).gunzip.map{it[1]}
            }
            else{
                ch_gtf = Channel.value(file(gtf))
            }
        }
        else if (gff) {
            if ( gff.endsWith(".gz")){
                ch_gff = GUNZIP_GTF([[:],gff]).gunzip.map{it[1]}
            }
            else{
                ch_gff = Channel.value(file(gff))
            }
            ch_gtf = GFFREAD ( ch_gff ).gtf
        }
    }

    CORRECT_INPUT(
        ch_fasta,
        ch_gtf
    )


    // ungzip adapter_fasta.gz file
    //if (adapter_fasta){
    //    if ( adapter_fasta.endsWith(".gz") ){
    //    ch_adapter_fasta = GUNZIP_ADAPTER_FASTA([[:],adapter_fasta]).gunzip.map{it[1]}
    //}else{
    //    ch_adapter_fasta = Channel.value(file(adapter_fasta))
   // }
    //}

    emit:
    ch_fasta = CORRECT_INPUT.out.out_fasta
    ch_gff = CORRECT_INPUT.out.out_gff
    //ch_adapter_fasta = ch_adapter_fasta
    
}
