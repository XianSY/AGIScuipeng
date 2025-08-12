include { FUSION_COMPUTE_WEIGHTS } from '../../../modules/local/fusion_compute_weights'
include { PLINK_FORMAT           } from '../../../modules/local/plink_format'
include { PHENOTYPE_DEAL         } from '../../../modules/local/phenotype_deal'
include { EXTRACT_GWAS_SUMMARY   } from '../../../modules/local/extract_gwas_summary'
include { MAKE_LD                } from '../../../modules/local/make_LD'
include { MERGE_WEIGHTS_INFO     } from '../../../modules/local/merge_weights_info'
include { PLINK_CHR_FORMAT       } from '../../../modules/local/plink_chr_format'
include { FUSION_ASSOC_TEST      } from '../../../modules/local/fusion_assoc_test'
include { SEPARATE_SNP_FOR_GENE  } from '../../../modules/local/separate_snp_for_gene'
include { GFF2GTF                } from '../../../modules/local/gff2gtf'

include { MAKE_GENE_INFO as MAKE_GENE_INFO_TWAS        } from '../../../modules/local/make_gene_info'



workflow FUSION_TWAS{

    main:
      
       ch_make_gene_info = [["id":"twas"],params.gff]
        if(params.gff.endsWith(".gff3") || params.gff.endsWith(".gff")){
            GFF2GTF(
                ch_make_gene_info
            )
            MAKE_GENE_INFO_TWAS(
            GFF2GTF.out.gtf
        )
        }
        else{
            MAKE_GENE_INFO_TWAS(
            ch_make_gene_info)
        }

        ch_extract_gwas = Channel.of([["id":"twas"],params.gwas_results,params.gwas_vcf])
   
        EXTRACT_GWAS_SUMMARY(ch_extract_gwas)        

 
       SEPARATE_SNP_FOR_GENE(
          [[id:"twas"],[params.eqtl_results]] 
       )
        PLINK_FORMAT(
            [[id:"twas"],[params.twas_compute_vcf]],
            SEPARATE_SNP_FOR_GENE.out.snp_dir 
        )

	//COMPUTE_ALLETE_RATE(
        //    [[id:"group"],[params.gwas_results_summary]]
        //)

        PLINK_CHR_FORMAT(
            [[id:"twas","chrome_num":params.chrome_num],[params.gwas_vcf]]
        )
        PHENOTYPE_DEAL(
           [[id:"twas"],[params.twas_compute_expression]]
        )
        	
	ch_compute_weights_parameter=Channel.of(["id":"twas"]).combine(PHENOTYPE_DEAL.out.expression)
        //ch_compute_weights_parameter.view()       
	FUSION_COMPUTE_WEIGHTS(
            ch_compute_weights_parameter,
            PLINK_FORMAT.out.plink_format
        )

	
	FUSION_COMPUTE_WEIGHTS.out.fusion_weights.collect().toList().set{ch_fusion_weight}
	 
        ch_fusion_weight.view()


        ch_gene_wgt_pos = Channel.of(["id":"twas"]).combine(FUSION_COMPUTE_WEIGHTS.out.fusion_weights)
	
        //ch_gene_wgt_pos.view()
	//ch_gene_wgt_pos = FUSION_COMPUTE_WEIGHTS.out.fusion_weights.map{
        //    it -> [it,params.gene_info]
        //}.groupTuple()

        MERGE_WEIGHTS_INFO(
            ch_gene_wgt_pos,
            MAKE_GENE_INFO_TWAS.out.gene_info
        )

	
        //把sumstats与染色体数及meta合并
        
        ch_summary_chr = Channel.of([["id":"twas","chrome_num":params.chrome_num],[params.gwas_results_summary]])	
        FUSION_ASSOC_TEST(
            PLINK_CHR_FORMAT.out.ldref,
       	    ch_summary_chr,
            MERGE_WEIGHTS_INFO.out.WEIGHT,
            MERGE_WEIGHTS_INFO.out.gene_wgt
        )
}
