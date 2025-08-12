include { ANNO_GFF3TOGENEPRED           } from '../../../modules/local/anno_gff3togenepred'
include { ANNO_RETRIEVE_SEQ             } from '../../../modules/local/anno_retrieve_seq_from_fasta'
include { ANNO_CONVERT2ANNOVAR          } from '../../../modules/local/anno_convert2annovar'
include { ANNO_TABLE                    } from '../../../modules/local/anno_table_annovar'



workflow ANNOTATION {
    
    take:
    ch_vcf_tbi_annotable
    ch_gff
    ch_fasta
    main:
    
    ch_anno_gff_meta=['id':'ANNO']
    ch_gff_gff3togenepred_withmeta=Channel.of(ch_anno_gff_meta).combine(ch_gff).collect().map{it}
    ch_gff_gff3togenepred_withoutmeta=ch_gff
    //ch_gff_gff3togenepred_withmeta=tuple(ch_anno_gff_meta, file(params.gff))
    //ch_gff_gff3togenepred_withoutmeta=file(params.gff)
    ANNO_GFF3TOGENEPRED(
        ch_gff_gff3togenepred_withmeta
    )

    //ch_fasta_gff3togenepred_withmeta=tuple(ch_anno_gff_meta, file(params.fasta))
    //ch_fasta_genepred_gff3togenepred_withmeta=ANNO_GFF3TOGENEPRED.out.genep.map{
    //    it-> tuple(it[0],it[1],ch_fasta_gff3togenepred_withmeta[1])}
    
    ch_fasta_gff3togenepred_withmeta=Channel.of(ch_anno_gff_meta).combine(ch_fasta).collect().map{it}
    ch_fasta_genepred_gff3togenepred_withmeta=ANNO_GFF3TOGENEPRED.out.genep.join(ch_fasta_gff3togenepred_withmeta)    
    ANNO_RETRIEVE_SEQ(
        ch_fasta_genepred_gff3togenepred_withmeta
    )
    
    ch_vcf_tbi_annotable=ch_vcf_tbi_annotable
    ANNO_CONVERT2ANNOVAR(
        ch_vcf_tbi_annotable
    )

   
    //ch_vcf_tbi_fa_fa2_gff_list_annotable=ch_vcf_tbi_annotable.map{
   //     it-> tuple(it[0],it[1],it[2],ch_fasta_gff3togenepred_withmeta[1])}.join(
   //                                 ANNO_RETRIEVE_SEQ.out.fa).map{
   //     it-> tuple(it[0],it[1],it[2],it[3],it[4],ch_gff_gff3togenepred_withmeta[1])}.join(
   //                                  ANNO_GFF3TOGENEPRED.out.genep).join(
   //                                     ANNO_CONVERT2ANNOVAR.out.anno
   //                                  )

   ch_vcf_tbi_fa_fa2_gff_list_annotable=ch_vcf_tbi_annotable.join(ch_fasta_gff3togenepred_withmeta)
                                        .join(ANNO_RETRIEVE_SEQ.out.fa)
                                        .join(ch_gff_gff3togenepred_withmeta)
                                        .join(ANNO_GFF3TOGENEPRED.out.genep)
                                        .join(ANNO_CONVERT2ANNOVAR.out.anno) 

    ANNO_TABLE(
        ch_vcf_tbi_fa_fa2_gff_list_annotable
    )
   
   emit:
   ch_annoevf_emmax_anno=ANNO_TABLE.out.annoevf.map{it->it[1]}

}
