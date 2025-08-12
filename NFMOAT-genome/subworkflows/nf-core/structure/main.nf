include { STURCTURE_VCF_ADDID           } from '../../../modules/local/structure_vcfaddid'
include { PLINK_INDEP                   } from '../../../modules/nf-core/plink/indep/main'
include { PLINK_EXTRACT                 } from '../../../modules/nf-core/plink/extract/main'
include { ADMIXTURE                     } from '../../../modules/nf-core/admixture/main'


workflow STRUCTURE {
    take:
    ch_vcf_tbi_structure_vcfaddid

    main:
    
    ch_vcf_tbi_structure_vcfaddid=ch_vcf_tbi_structure_vcfaddid
   
    
    STURCTURE_VCF_ADDID(
        ch_vcf_tbi_structure_vcfaddid
    )
    
    ch_vcf_tbi_plink_indep=STURCTURE_VCF_ADDID.out.vcf

    PLINK_INDEP(
       ch_vcf_tbi_plink_indep
    )
    plink_indep_prunein = PLINK_INDEP.out.prunein
    plink_indep_prunein.view()
    ch_vcf_tbi_prunein_plink_extract=ch_vcf_tbi_plink_indep.join(plink_indep_prunein)
    
    PLINK_EXTRACT(
        ch_vcf_tbi_prunein_plink_extract
    )
    
    ch_bed_bim_fam_k_admixture=PLINK_EXTRACT.out.bed.join(PLINK_EXTRACT.out.bim).join(PLINK_EXTRACT.out.fam).combine(Channel.from(2,3,4,5,6,7,8,9,10)).map{
        it -> tuple(['id':"popu_snp_${it[4]}"],it[1],it[2],it[3],it[4])
    }
    //ch_bed_bim_fam_k_admixture.view()
    ADMIXTURE(
        ch_bed_bim_fam_k_admixture
    )
            //    R_STRUCTURE_PLOT()
    
   emit:
   ch_plink_indep_prunein = plink_indep_prunein
}
