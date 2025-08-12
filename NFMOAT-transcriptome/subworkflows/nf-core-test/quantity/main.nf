include { FEATURECOUNTS                      } from '../../../modules/local/featurecounts'
include { FILTER_GENE_EXPRESSION             } from '../../../modules/local/filter_gene_expression'
include { TPM_PCA                            } from '../../../modules/local/tpm_pca'
include { QUALIMAP_RNASEQ                    } from '../../../modules/nf-core/qualimap/rnaseq/main'

workflow QUANTIFY{

    take:
        gtf
        ch_genome_bam_total
        ch_genome_bam
	biotype        

    main:

        ch_qualimap_multiqc           = Channel.empty()
        //
        QUALIMAP_RNASEQ (
                ch_genome_bam,
                gtf
            )
            ch_qualimap_multiqc = QUALIMAP_RNASEQ.out.results


	//ch_genome_bam.view()
        //
        // MODULE: Feature biotype QC using featureCounts
        //
        def meta_total=["group":"counts"]
        ch_gtf_featurecounts = gtf.map{
            it -> tuple(meta_total,it)
        }.groupTuple()
        ch_featurecounts_multiqc = Channel.empty()
        if (!params.skip_alignment && !params.skip_qc && !params.skip_biotype_qc && biotype) {

            gtf
            .map { WorkflowRnaseq.biotypeInGtf(it, biotype, log) }
            .set { biotype_in_gtf }

            // Prevent any samples from running if GTF file doesn't have a valid biotype
            ch_genome_bam
                .combine(gtf)
                .combine(biotype_in_gtf)
                .filter { it[-1] }
                .map { it[0..<it.size()-1] }
                .set { ch_featurecounts }
            FEATURECOUNTS (
            ch_gtf_featurecounts,  
            ch_genome_bam_total
            ) 
            FILTER_GENE_EXPRESSION(
                FEATURECOUNTS.out.counts.map{it -> tuple(it[0],it[1])}.groupTuple()
            )
            ch_deseq_featurecounts = FILTER_GENE_EXPRESSION.out.filter_gene_expression.map{
              it -> tuple(it[0],it[1],params.input_compare,params.input_group)
              }.groupTuple()
            ch_tpm_pca = FEATURECOUNTS.out.featCounts.map{it -> tuple(it[0],it[1],params.input_group)}.groupTuple()
 	   
        } 
            ch_tpm_pca.view()
   	    TPM_PCA(
                ch_tpm_pca
            )
           ch_tpm = TPM_PCA.out.deseq.map{it -> tuple(it[0],it[1],it[2])}.groupTuple() 

    emit:

        ch_deseq_featurecounts
        ch_tpm_pca
        ch_featurecounts
        ch_tpm


}
