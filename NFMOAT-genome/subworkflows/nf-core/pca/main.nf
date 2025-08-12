include { PCA_VCF_ADDID                 } from '../../../modules/local/pca_vcfaddid'
include { VCFTOOLS_PLINK                } from '../../../modules/local/vcftools_plink'
include { PLINK_MAKEBED                 } from '../../../modules/local/plink_makebed'
include { GCTA_MAKEGRM                  } from '../../../modules/local/gcta_makegrm'
include { GCTA_PCA                      } from '../../../modules/local/gcta_pca'
include { GCTA_PCA_PLOT                 } from '../../../modules/local/gcta_pca_plot'


workflow PCA {
    take:
    ch_vcf_tbi_pca_vcfaddid
    ch_plink_indep_prunein
    
    main:

    PCA_VCF_ADDID(
        ch_vcf_tbi_pca_vcfaddid
    )
    
    ch_vcftools_plink=PCA_VCF_ADDID.out.vcf.join(ch_plink_indep_prunein)
   
    VCFTOOLS_PLINK(
        ch_vcftools_plink
    )

    PLINK_MAKEBED(
        VCFTOOLS_PLINK.out.plink
    )

    ch_bed_bim_fam_gcta_makegrm=VCFTOOLS_PLINK.out.plink.join(PLINK_MAKEBED.out.bed)
    GCTA_MAKEGRM(
        ch_bed_bim_fam_gcta_makegrm
    )

    GCTA_PCA(
        GCTA_MAKEGRM.out.grm
    )
    ch_eigenvec = GCTA_PCA.out.vec.map{
        row -> tuple(row[0],row[2])
    }
    GCTA_PCA_PLOT(
        ch_eigenvec
    )

}
