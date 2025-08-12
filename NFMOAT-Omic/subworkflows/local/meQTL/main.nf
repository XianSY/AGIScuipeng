include { BETA_MATRIX_TO_MEQTL    } from '../../../modules/local/beta_matrix_to_meqtl'
include { VCF_TO_NUMBER          } from '../../../modules/local/vcf_to_num'
include { VCF_TO_EQTL  as VCF_TO_MEQTL          } from '../../../modules/local/vcf_to_eqtl'
include { EQTL_PROCESS  as MEQTL_PROCESS } from '../../../modules/local/eqtl'
include { VCF_TO_INDENTIFY  as VCF_TO_MEQTL_INDENTIFY      } from '../../../modules/local/vcf_to_indentify'
include { VCF_TO_NUM_EQTL as VCF_TO_NUM_MEQTL        } from '../../../modules/local/vcf_to_num_eqtl'
include { VCF_TO_BIALLELIC as VCF_TO_MEQTL_BIALLELIC      } from '../../../modules/local/vcf_to_biallelic'


workflow MEQTL{

    main:

        ch_beta_matrix = [["id":"meQTL"],params.beta_matrix]

        BETA_MATRIX_TO_MEQTL(
            ch_beta_matrix
        )        
 
        //snp file for genome variant calling 

        ch_vcf_to_indetify = [["id":"meQTL"],params.meQTL_vcf]

        //VCF_TO_NUMBER(
        //    ch_vcf_to_num
        //)
        
        //确定vcf的格式
        VCF_TO_MEQTL_INDENTIFY(
            ch_vcf_to_indetify
        )
        
        //去除复等位基因
        VCF_TO_MEQTL_BIALLELIC(
            VCF_TO_MEQTL_INDENTIFY.out.gt_vcf
        )

        VCF_TO_NUM_MEQTL(
            VCF_TO_MEQTL_BIALLELIC.out.biallelic_vcf
        ) 

        //VCF_TO_MEQTL(
        //    ch_vcf_to_num
        //)

        ch_snp_eqtl = VCF_TO_NUM_MEQTL.out.vcf_num.flatten().collect()

        ch_beta_matrix_eqtl = BETA_MATRIX_TO_MEQTL.out.smp_meQTL

        ch_snp_eqtl.view()
        ch_beta_matrix_eqtl.view()

        MEQTL_PROCESS(
          ch_snp_eqtl,
          ch_beta_matrix_eqtl
        )
}
