include { BEDGRAPH2BOD as BEDGRAPH2BOD_CHG     } from '../../../modules/local/bedgraph2bod.nf'
include { BEDGRAPH2BOD as BEDGRAPH2BOD_CPG     } from '../../../modules/local/bedgraph2bod.nf'
include { BEDGRAPH2BOD as BEDGRAPH2BOD_CHH     } from '../../../modules/local/bedgraph2bod.nf'
include { BEDTRANFORMBOD                       } from '../../../modules/local/bedtranfrombod.nf'
include { OSCAPCA as OSCAPCA_CHG               } from '../../../modules/local/osca_pca.nf'
include { OSCAPCA as OSCAPCA_CPG               } from '../../../modules/local/osca_pca.nf'
include { OSCAPCA as OSCAPCA_CHH               } from '../../../modules/local/osca_pca.nf'
include { OSCAEWAS as OSCAEWASPCA_CHG       } from '../../../modules/local/osca_ewas.nf'
include { OSCAEWAS as OSCAEWASPCA_CPG       } from '../../../modules/local/osca_ewas.nf'
include { OSCAEWAS as OSCAEWASPCA_CHH       } from '../../../modules/local/osca_ewas.nf'
include { PLOT_HISTORGAM as PLOT_HISTORGAM_CHG                  } from '../../../modules/local/plot_histogram'
include { PLOT_HISTORGAM as PLOT_HISTORGAM_CPG                  } from '../../../modules/local/plot_histogram'
include { PLOT_HISTORGAM as PLOT_HISTORGAM_CHH                  } from '../../../modules/local/plot_histogram'
include { BETA2MENOYPTE  as BETE2MENOYPTE_CHG                } from '../../../modules/local/beta2menoypte'
include { BETA2MENOYPTE  as BETE2MENOYPTE_CPG                } from '../../../modules/local/beta2menoypte'
include { BETA2MENOYPTE  as BETE2MENOYPTE_CHH                } from '../../../modules/local/beta2menoypte'
include { MENOYPTE2VCF  as MENOYPTE2VCF_CHG                } from '../../../modules/local/menoypte2vcf'
include { MENOYPTE2VCF  as MENOYPTE2VCF_CPG                } from '../../../modules/local/menoypte2vcf'
include { MENOYPTE2VCF  as MENOYPTE2VCF_CHH                } from '../../../modules/local/menoypte2vcf'
include { SNPEFF_ANN_VCF as SNPEFF_ANN_VCF_CHG                  } from '../../../modules/local/snpEff_ann_vcf'
include { SNPEFF_ANN_VCF as SNPEFF_ANN_VCF_CPG                  } from '../../../modules/local/snpEff_ann_vcf'
include { SNPEFF_ANN_VCF as SNPEFF_ANN_VCF_CHH                  } from '../../../modules/local/snpEff_ann_vcf'
include { PLINK_MAKEBED  as PLINK_MAKEBED_CPG                   } from '../../../modules/local/plink_make_bed'
include { PLINK_MAKEBED  as PLINK_MAKEBED_CHG                   } from '../../../modules/local/plink_make_bed'
include { PLINK_MAKEBED  as PLINK_MAKEBED_CHH                   } from '../../../modules/local/plink_make_bed'
include { PLINK_PRE_MAKEBED as PLINK_PRE_MAKEBED_CPG            } from '../../../modules/local/plink_pre_make_bed'
include { PLINK_PRE_MAKEBED as PLINK_PRE_MAKEBED_CHG            } from '../../../modules/local/plink_pre_make_bed'
include { PLINK_PRE_MAKEBED as PLINK_PRE_MAKEBED_CHH            } from '../../../modules/local/plink_pre_make_bed'
include { EWAS_EMMAX_KINSHIP as EWAS_EMMAX_KINSHIP_CPG          } from '../../../modules/local/Ewas_emmax_kin'
include { EWAS_EMMAX_KINSHIP as EWAS_EMMAX_KINSHIP_CHG          } from '../../../modules/local/Ewas_emmax_kin'
include { EWAS_EMMAX_KINSHIP as EWAS_EMMAX_KINSHIP_CHH          } from '../../../modules/local/Ewas_emmax_kin'
include { STURCTURE_VCF_ADDID as STURCTURE_VCF_ADDID_CPG        } from  '../../../modules/local/structure_vcfadd'
include { STURCTURE_VCF_ADDID as STURCTURE_VCF_ADDID_CHG        } from  '../../../modules/local/structure_vcfadd'
include { STURCTURE_VCF_ADDID as STURCTURE_VCF_ADDID_CHH        } from  '../../../modules/local/structure_vcfadd'
include { GETCHRLENGTH                                          } from '../../../modules/local/getchrlength'
include { VCF_ADD_HEADER as VCF_ADD_HEADER_CHG                  } from '../../../modules/local/vcf_add_header'
include { VCF_ADD_HEADER as VCF_ADD_HEADER_CPG                  } from '../../../modules/local/vcf_add_header'
include { VCF_ADD_HEADER as VCF_ADD_HEADER_CHH                  } from '../../../modules/local/vcf_add_header'
include { PLINK_INDEP as EMMAX_PLINK_INDEP                      } from '../../../modules/nf-core/plink/indep/main'
include { PLINK_EXTRACT as EMMAX_PLINK_EXTRACT_CPG              } from '../../../modules/nf-core/plink/extract/main'
include { PLINK_INDEP   as EMMAX_PLINK_INDEP_CPG                } from '../../../modules/nf-core/plink/indep/main'
include { ADMIXTURE as EMMAX_ADMIXTURE_CPG                      } from '../../../modules/nf-core/admixture/main'
include { GWAS_EMMAX_COV  as EWAS_EMMAX_COV_CPG                 } from '../../../modules/local/gwas_emmax_cov'

workflow OSCAEWAS {

    take:
    bedgraph  
    chg_bedgraph  //这个是要转换成ewas的beta—value文件
    cpg_bedgraph
    chh_bedgraph
    reference_index_fai
    main:
  
    //PLOT_HISTORGAM_CHG(
    //    chg_bedgraph,
    //    'chg'
   // )
   // PLOT_HISTORGAM_CPG(
   //     cpg_bedgraph,
   //     'cpg'
   // )
   // PLOT_HISTORGAM_CHH(
   //     chh_bedgraph,
   //     'chh'
   // )    

    

    BEDGRAPH2BOD_CHG(
       chg_bedgraph
    )
    BEDGRAPH2BOD_CPG(
        cpg_bedgraph
    )
    BEDGRAPH2BOD_CHH(
        chh_bedgraph
    )
   
   BEDGRAPH2BOD_CHG.out.bod.view()
 
    BETE2MENOYPTE_CHG(
        chg_bedgraph
    )
    //println(BETE2MENOYPTE_CHG.out.smp_menoypte)
    BETE2MENOYPTE_CHG.out.smp_menoypte.view()
    
    BETE2MENOYPTE_CPG(
        cpg_bedgraph
    )
    BETE2MENOYPTE_CHH(
        chh_bedgraph
    )
    //ch_smp_menoypte = BETE2MENOYPTE_CHG.out.menoypte
    //println(111111111111111)
    //println(ch_smp_menoypte)
    //BETE2MENOYPTE_CHG.out.smp_menoypte.view()   

    MENOYPTE2VCF_CHG(
        BETE2MENOYPTE_CHG.out.smp_menoypte
    )
    MENOYPTE2VCF_CPG(
        BETE2MENOYPTE_CPG.out.smp_menoypte
    )
    MENOYPTE2VCF_CHH(
        BETE2MENOYPTE_CHH.out.smp_menoypte
    )
    
    //MENOYPTE2VCF_CHG.out.SMP_vcf.view()
   // println(MENOYPTE2VCF_CHG.out.SMP_vcf)       
    
 
    chg_menotype = MENOYPTE2VCF_CHG.out.SMP_vcf.map{row -> tuple(row[0],row[1],params.fasta,params.gff)}.groupTuple().view()
    cpg_menotype = MENOYPTE2VCF_CPG.out.SMP_vcf.map{row -> tuple(row[0],row[1],params.fasta,params.gff)}.groupTuple().view()
    chh_menotype = MENOYPTE2VCF_CHH.out.SMP_vcf.map{row -> tuple(row[0],row[1],params.fasta,params.gff)}.groupTuple().view()    

   //tuple([MENOYPTE2VCF_CHG.out.SMP_vcf,[params.fasta],[params.gff]])
   // cpg_menotype = tuple([MENOYPTE2VCF_CPG.out.SMP_vcf,[params.fasta],[params.gff]])
   // chh_menotype = tuple([MENOYPTE2VCF_CHH.out.SMP_vcf,[params.fasta],[params.gff]])
    println(chg_menotype)
    SNPEFF_ANN_VCF_CHG(
        chg_menotype
    )
    SNPEFF_ANN_VCF_CPG(
        cpg_menotype
    )
    SNPEFF_ANN_VCF_CHH(
        chh_menotype
    )
    GETCHRLENGTH(
        reference_index_fai
    )
    
    VCF_ADD_HEADER_CHG(
        MENOYPTE2VCF_CHG.out.SMP_vcf,
        GETCHRLENGTH.out.contig
    )
    VCF_ADD_HEADER_CPG(
        MENOYPTE2VCF_CPG.out.SMP_vcf,
        GETCHRLENGTH.out.contig
    )
    VCF_ADD_HEADER_CHH(
        MENOYPTE2VCF_CHH.out.SMP_vcf,
        GETCHRLENGTH.out.contig
    )

    STURCTURE_VCF_ADDID_CHG(
        VCF_ADD_HEADER_CHG.out.SMP_vcf
    )
    STURCTURE_VCF_ADDID_CHH(
        VCF_ADD_HEADER_CHH.out.SMP_vcf
    )
    STURCTURE_VCF_ADDID_CPG(
        VCF_ADD_HEADER_CPG.out.SMP_vcf
    )
 
    EMMAX_PLINK_INDEP_CPG(
        STURCTURE_VCF_ADDID_CPG.out.vcf
    )

    ch_vcf_tbi_prunein_emmax_plink_extract_cpg=STURCTURE_VCF_ADDID_CPG.out.vcf.join(EMMAX_PLINK_INDEP_CPG.out.prunein)
    EMMAX_PLINK_EXTRACT_CPG(
        ch_vcf_tbi_prunein_emmax_plink_extract_cpg
    )
   
    ch_bed_bim_fam_EMMAX_k_admixture_cpg=EMMAX_PLINK_EXTRACT_CPG.out.bed.join(EMMAX_PLINK_EXTRACT_CPG.out.bim).join(EMMAX_PLINK_EXTRACT_CPG.out.fam).combine(Channel.from(2,3,4,5,6,7,8,9,10)).map{
        it -> tuple(['id':"popu_snp_${it[4]}","group":"${it[0].id}"],it[1],it[2],it[3],it[4])
    }
    //ch_bed_bim_fam_EMMAX_k_admixture.view()
    EMMAX_ADMIXTURE_CPG(
        ch_bed_bim_fam_EMMAX_k_admixture_cpg
    )
 
    PLINK_PRE_MAKEBED_CPG(
        VCF_ADD_HEADER_CPG.out.SMP_vcf
    )
    PLINK_PRE_MAKEBED_CHG(
        VCF_ADD_HEADER_CHG.out.SMP_vcf
    )
    PLINK_PRE_MAKEBED_CHH(
        VCF_ADD_HEADER_CHH.out.SMP_vcf
    )
    
    PLINK_MAKEBED_CPG(
        PLINK_PRE_MAKEBED_CPG.out.bed.map{
            row -> tuple(row[0],row[1],row[2],row[3])
        }
    )
    PLINK_MAKEBED_CHG(
        PLINK_PRE_MAKEBED_CHG.out.bed.map{
            row -> tuple(row[0],row[1],row[2],row[3])
        }
    )
    PLINK_MAKEBED_CHH(
        PLINK_PRE_MAKEBED_CHH.out.bed.map{
            row -> tuple(row[0],row[1],row[2],row[3])
        }
    )

    EWAS_EMMAX_KINSHIP_CPG(
        PLINK_MAKEBED_CPG.out.emmaxfile
    )
    EWAS_EMMAX_KINSHIP_CHG(
        PLINK_MAKEBED_CHG.out.emmaxfile
    )
    EWAS_EMMAX_KINSHIP_CHH(
        PLINK_MAKEBED_CHH.out.emmaxfile
    )

    ch_qfile_tpedfam_emmax_cov_cpg=EMMAX_ADMIXTURE_CPG.out.ancestry_fractions.join(EMMAX_ADMIXTURE_CPG.out.cvlog).map{
              it-> tuple(['id':"${it[0].group}"],it[1],it[2])}.groupTuple().join(
                PLINK_MAKEBED_CPG.out.emmaxfile.map{it-> tuple(['id':"${it[0].id}"],it[2])})
    EWAS_EMMAX_COV_CPG(
        ch_qfile_tpedfam_emmax_cov_cpg
    )

    ch_bedbod = BEDGRAPH2BOD_CHG.out.bod.join(BEDGRAPH2BOD_CPG.out.bod).join(BEDGRAPH2BOD_CHH.out.bod).view()
    BEDTRANFORMBOD(
        ch_bedbod
    )
    
    ch_oscapca_chg_bod = BEDTRANFORMBOD.out.CHG_BOD
       
    OSCAPCA_CHG(
        BEDTRANFORMBOD.out.CHG_BOD
    )
   OSCAPCA_CPG(
        BEDTRANFORMBOD.out.CPG_BOD
    )
    OSCAPCA_CHH(
        BEDTRANFORMBOD.out.CHH_BOD
    )
    
    //OSCAEWASPCA_CHG(ch_chg_ewas)
    //OSCAEWASPCA_CPG(ch_cpg_ewas)
    //OSCAEWASPCA_CHH(ch_chh_ewas)


    
}
