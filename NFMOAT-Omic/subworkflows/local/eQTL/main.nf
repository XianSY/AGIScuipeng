include { MAKE_GENE_INFO         } from '../../../modules/local/make_gene_info'
include { GENE_EXPRESSION_FILTER } from '../../../modules/local/gene_expression_filter'
include { VCF_TO_NUMBER          } from '../../../modules/local/vcf_to_num'
include { VCF_TO_EQTL            } from '../../../modules/local/vcf_to_eqtl'
include { EQTL_PROCESS           } from '../../../modules/local/eqtl'
include { VCF_TO_INDENTIFY       } from '../../../modules/local/vcf_to_indentify'
include { VCF_TO_NUM_EQTL        } from '../../../modules/local/vcf_to_num_eqtl'
include { VCF_TO_BIALLELIC       } from '../../../modules/local/vcf_to_biallelic'
include { GFF2GTF                } from '../../../modules/local/gff2gtf'
include { EQTL_RESULTS_FILTER    } from '../../../modules/local/eqtl_results_filter'


workflow EQTL{

    main:
        ch_make_gene_info = [["id":"eQTL"],params.gff]
        if(params.gff.endsWith(".gff3") || params.gff.endsWith(".gff")){
            GFF2GTF(
                ch_make_gene_info
            )
            MAKE_GENE_INFO(
            GFF2GTF.out.gtf
        )
        }
        else{
            MAKE_GENE_INFO(
            ch_make_gene_info)
        }   
     
        ch_gene_expression_filter = MAKE_GENE_INFO.out.gene_info.combine(Channel.fromPath(params.input_expression_eQTL))
        GENE_EXPRESSION_FILTER(
            ch_gene_expression_filter
        )
	ch_vcf_to_indetify = [["id":"eQTL"],params.input_vcf]
        //VCF_TO_NUMBER(
        //    ch_vcf_to_num
        //)
        //确定vcf的格式
        VCF_TO_INDENTIFY(
            ch_vcf_to_indetify
        )


        VCF_TO_BIALLELIC(
            VCF_TO_INDENTIFY.out.gt_vcf
        )
        
        VCF_TO_NUM_EQTL(
            VCF_TO_BIALLELIC.out.biallelic_vcf
        )         


	ch_snp_eqtl = VCF_TO_NUM_EQTL.out.vcf_num
        ch_expression_eqtl = GENE_EXPRESSION_FILTER.out.gene_expression_eQTL
	EQTL_PROCESS(
	  ch_snp_eqtl,
	  ch_expression_eqtl,
	  VCF_TO_INDENTIFY.out.eigenvec
	)
	EQTL_RESULTS_FILTER(
            ch_snp_eqtl,
            EQTL_PROCESS.out.init_eqtl,
	    ch_expression_eqtl

        )
}
