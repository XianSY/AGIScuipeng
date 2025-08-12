include { MAKE_GENE_INFO as MAKE_GENE_INFO_EQTM         } from '../../../modules/local/make_gene_info'
include { GENE_EXPRESSION_FILTER as GENE_EXPRESSION_FILTER_EQTM } from '../../../modules/local/gene_expression_filter'
include { VCF_TO_NUMBER          } from '../../../modules/local/vcf_to_num'
include { VCF_TO_EQTM            } from '../../../modules/local/smp_vcf_to_eqtm'
include { EQTL_PROCESS as EQTM_PROCESS          } from '../../../modules/local/eqtl'
include { VCF_TO_INDENTIFY as VCF_TO_MENOTYPE_INDENTIFY       } from '../../../modules/local/vcf_to_indentify'
include { VCF_TO_NUM_EQTL  as VCF_TO_NUM_EQTM      } from '../../../modules/local/vcf_to_num_eqtl'
include { VCF_TO_BIALLELIC as VCF_TO_MENOTYPE_BIALLELIC      } from '../../../modules/local/vcf_to_biallelic'
include { GFF2GTF                } from '../../../modules/local/gff2gtf'


workflow EQTM{

    main:
        ch_make_gene_info = [["id":"eQTM"],params.gff]

        if(params.gff.endsWith(".gff3") || params.gff.endsWith(".gff")){
            GFF2GTF(
                ch_make_gene_info
            )
            MAKE_GENE_INFO_EQTM(
            GFF2GTF.out.gtf
        )
        }
        else{
            MAKE_GENE_INFO_EQTM(
            ch_make_gene_info)
        }   

        ch_gene_expression_filter = MAKE_GENE_INFO_EQTM.out.gene_info.combine(Channel.fromPath(params.input_expression_eQTM))
        GENE_EXPRESSION_FILTER_EQTM(
            ch_gene_expression_filter
        )
	
	ch_vcf_to_indetify = [["id":"eQTM"],params.smp_vcf]
       
        VCF_TO_MENOTYPE_INDENTIFY(
            ch_vcf_to_indetify
        )

        VCF_TO_MENOTYPE_BIALLELIC(
            VCF_TO_MENOTYPE_INDENTIFY.out.gt_vcf
        )
        
	VCF_TO_NUM_EQTM(
            VCF_TO_MENOTYPE_BIALLELIC.out.biallelic_vcf
        )

        ch_smp_eqtl = VCF_TO_NUM_EQTM.out.vcf_num
        ch_expression_eqtm = GENE_EXPRESSION_FILTER_EQTM.out.gene_expression_eQTL
        EQTM_PROCESS(
          ch_smp_eqtl,
          ch_expression_eqtm
        )        

}
