include { DESEQ2_COUNTS                      } from '../../../modules/local/deseq2_featurecounts'
include { EDGER_COUNTS                       } from '../../../modules/local/edgeR_featurecount'
include { GET_ENRICEMENG_GENE                } from '../../../modules/local/get_gene'
include { GO_KEGG  as GO_KEGG_DESEQ2         } from '../../../modules/local/go_kegg'
include { GO_KEGG  as GO_KEGG_EDGER          } from '../../../modules/local/go_kegg'

workflow DIFFERENCE_EXPRESSION{
    take:
        ch_deseq_featurecounts
        ch_tpm_pca
    main:
	//println(ch_deseq_featurecounts)
        //println(ch_tpm_pca)
	DESEQ2_COUNTS(
                ch_deseq_featurecounts
            )
        
        EDGER_COUNTS(
                ch_deseq_featurecounts
            )
	if(params.species != ''){
	       GET_ENRICEMENG_GENE(DESEQ2_COUNTS.out.deseq_results.flatten())
               GO_KEGG_DESEQ2(GET_ENRICEMENG_GENE.out.get_gene.map{row -> tuple(['id':'enrichment'],row,params.species)})

	}
}
