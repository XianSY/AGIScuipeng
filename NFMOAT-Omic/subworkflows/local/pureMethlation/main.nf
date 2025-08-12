include { FAMDB        } from '../../../modules/local/famdb'
include { REPEATMASKER } from '../../../modules/local/repeatmasker'  
include { MAKE_BED as MAKE_BED_METHYLATION     } from '../../../modules/local/make_bed'
include { MAKE_BED as MAKE_BED_SNP             } from '../../../modules/local/make_bed'
include { MAKE_BED as MAKE_BED_INDEL           } from '../../../modules/local/make_bed'
include { MAKE_BED_REPEATMASKER                } from '../../../modules/local/make_bed_repeatmaskers'
include{ BEDTOOLS_INTER as BEDTOOLS_INTER_METHYLATION_SNP          } from '../../../modules/local/bedtools_inter'
include{ BEDTOOLS_INTER as BEDTOOLS_INTER_METHYLATION_INDEL        } from '../../../modules/local/bedtools_inter'
include{ BEDTOOLS_INTER as BEDTOOLS_INTER_METHYLATION_REPEATMASKER } from '../../../modules/local/bedtools_inter'


workflow PUREMETHLATION{
    main:
        ch_famdb_input = Channel.of([["id":"TE_anno"]])
        ch_famdb_input.view()
        //denovo annotation TE
        FAMDB(
            ch_famdb_input
        )

        ch_repeatmasker = FAMDB.out.ch_famdb.map{ it -> tuple(it[0],it[2])}.combine(Channel.fromPath(params.fasta))

        REPEATMASKER(
            ch_repeatmasker
        )

        MAKE_BED_REPEATMASKER(REPEATMASKER.out.fasta)
	ch_to_bed = Channel.of([["id":"PURE"],file(params.methylation),file(params.snp),file(params.InDel)])
                    .combine(REPEATMASKER.out.fasta)
                    .map{it -> tuple(it[0],it[1],it[2],it[3],it[5])}
        //ch_to_bed.view()
        
        ch_methylation_to_bed = Channel.of([["id":"PURE"],file(params.methylation)])
        ch_snp_to_bed = Channel.of([["id":"PURE"],file(params.snp)])
        ch_indel_to_bed = Channel.of([["id":"PURE"],file(params.InDel)])
        
        MAKE_BED_METHYLATION(ch_methylation_to_bed)
        MAKE_BED_SNP(ch_snp_to_bed)
        MAKE_BED_INDEL(ch_indel_to_bed)     
        ch_methylation_snp = MAKE_BED_METHYLATION.out.bed.join(MAKE_BED_SNP.out.bed)
        BEDTOOLS_INTER_METHYLATION_SNP(
	  ch_methylation_snp
        )   
        
        ch_methylation_indel = BEDTOOLS_INTER_METHYLATION_SNP.out.pure_methy_file.join(MAKE_BED_INDEL.out.bed)
        ch_methylation_indel.view()
        BEDTOOLS_INTER_METHYLATION_INDEL(
            ch_methylation_indel
        )        
	
        ch_methylation_repeatmasker = BEDTOOLS_INTER_METHYLATION_INDEL.out.pure_methy_file
                                      .combine(MAKE_BED_REPEATMASKER.out.bed).flatten().collect()
                                      .map{it -> tuple(it[0],it[1],it[3])}.view()
        //ch_methylation_repeatmasker.view()
        BEDTOOLS_INTER_METHYLATION_REPEATMASKER(ch_methylation_repeatmasker)
                
        
}
