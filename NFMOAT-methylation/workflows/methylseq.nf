/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap; fromSamplesheet } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowMethylseq.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOWS: Consisting of a mix of local and nf-core/modules
//
include { PREPARE_GENOME } from '../subworkflows/local/prepare_genome'
include { FASTQC_WORK    } from '../subworkflows/local/fastqc/fastqc'
// Aligner: bismark or bismark_hisat
if( params.aligner =~ /bismark/ ){
    include { BISMARK_SUB     } from '../subworkflows/local/align/bismark'
    include { BISMARK_DMR     } from '../subworkflows/local/dmr/bismark_dmr'
    include { OSCAEWAS        } from '../subworkflows/local/osca_ewas/oscaewas'
}
// Aligner: bwameth
else if ( params.aligner == 'bwameth' ){
    include { BWAMETH } from '../subworkflows/local/bwameth'
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { QUALIMAP_BAMQC              } from '../modules/nf-core/qualimap/bamqc/main'
include { PRESEQ_LCEXTRAP             } from '../modules/nf-core/preseq/lcextrap/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow METHYLSEQ{

    ch_versions = Channel.empty()

    processList = params.processList.split(',')
    start =  processList[0]
    //println(processList)
    //println(params.input)
    // SUBWORKFLOW: Prepare any required reference genome indices
    //

    PREPARE_GENOME() //为参考基因组建立索引，如果有索引就不建立
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)
    samtools_kary = PREPARE_GENOME.out.samtools_kary.map{row -> tuple([id:"group"],row[1],row[2])}.groupTuple()
    //ch_versions.view() 
    //
    // Create input channel from input file provided through params.input
    //
    if(start == "A"){
        Channel
        .fromSamplesheet("input")
        .map {
            meta, fastq_1, fastq_2 ->
            if (!fastq_2) {
                return [ meta + [ single_end:true ], [ fastq_1 ] ]
            } else {
                return [ meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
            }
        }
        .groupTuple()
        .map {
            meta, fastq ->
            def meta_clone = meta.clone()
            parts = meta_clone.id.split('_')
            meta_clone.id = parts.length > 1 ? parts[0..-2].join('_') : meta_clone.id
            [ meta_clone, fastq ]
        }
        .groupTuple(by: [0])
        .branch {
            meta, fastq ->
            single: fastq.size() == 1
            return [ meta, fastq.flatten() ]
            multiple: fastq.size() > 1
            return [ meta, fastq.flatten() ]
        }
        .set { ch_fastq }
     FASTQC_WORK()
     //FASTQC_WORK.out.reads.view() 
    }
    if(start == "B"){
            /*
        * SUBWORKFLOW: Align reads, deduplicate and extract methylation with Bismark
        */

        // Aligner: bismark or bismark_hisat
        reads = Channel.fromPath(params.input_align).splitCsv(header:true).map{
		row -> if(!row.fastq_2){
                         return [[id:row.sample]+[single_end:true],[row.fastq_1]]
               }else{
                  return [ [id:row.sample]+[single_end:false],[row.fastq_1,row.fastq_2]]
               }
        }
        .map {
            meta, fastq ->
            def meta_clone = meta.clone()
            parts = meta_clone.id.split('_')
            meta_clone.id = parts.length > 1 ? parts[0..-2].join('_') : meta_clone.id
            [ meta_clone, fastq.flatten() ]
        }
        //reads = Channel.fromSamplesheet("input_align")
        //.map {
        //    meta, fastq_1, fastq_2 ->
        //    if (!fastq_2) {
        //        return [ meta + [ single_end:true ], [ fastq_1 ] ]
        //    } else {
        //        return [ meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
        //     }
        // }
        //.groupTuple()
        //.map {
        //    meta, fastq ->
        //    def meta_clone = meta.clone()
        //    parts = meta_clone.id.split('_')
        //    meta_clone.id = parts.length > 1 ? parts[0..-2].join('_') : meta_clone.id
        //    [ meta_clone, fastq.flatten() ]
        //}

 	//reads = ch_fastq.multiple.view()
  
        if( params.aligner =~ /bismark/ ){

            /*
            * Run Bismark alignment + downstream processing
            */

            BISMARK_SUB (
                start,
                reads,
                PREPARE_GENOME.out.bismark_index,
                params.skip_deduplication || params.rrbs,
                params.cytosine_report || params.nomeseq,
                summary_params,
                ch_multiqc_custom_methods_description,
                ch_multiqc_config,
                ch_multiqc_custom_config,
                ch_multiqc_logo
            )
                  
           // ch_versions = ch_versions.mix(BISMARK_SUB.out.versions.unique{ it.baseName })
		
   }
        }else if('B' in processList){
        BISMARK_SUB(
                start,
                FASTQC_WORK.out.reads,
                PREPARE_GENOME.out.bismark_index,
                params.skip_deduplication || params.rrbs,
                params.cytosine_report || params.nomeseq,
                summary_params,
                ch_multiqc_custom_methods_description,
                ch_multiqc_config,
                ch_multiqc_custom_config,
                ch_multiqc_logo,
                )

            //ch_versions = ch_versions.mix(BISMARK.out.versions.unique{ it.baseName })
            //ch_bam = BISMARK.out.bam
            //ch_dedup = BISMARK.out.dedup
            // ch_aligner_mqc = BISMARK.out.mqc
    }
    if(start =="C"){ //进行DMR和EWAS
	ch_methylation_chg = Channel.fromPath(params.input_chg)
        .splitCsv( header:true,sep:",")
        .map{
            row -> [[id:row.sample],row.deduplicated,row.report]
        }
        .groupTuple()

        ch_methylation_cpg = Channel.fromPath(params.input_cpg)
        .splitCsv( header:true,sep:",")
        .map{
            row -> [[id:row.sample],row.deduplicated,row.report]
        }
        .groupTuple()

        ch_methylation_chh = Channel.fromPath(params.input_chh)
        .splitCsv( header:true,sep:",")
        .map{
            row -> [[id:row.sample],row.deduplicated,row.report]
        }
        .groupTuple()

       BISMARK_DMR(
   	 ch_methylation_chg,
	 ch_methylation_cpg,
     	 ch_methylation_chh
	 )

	}else if('C' in processList){
            BISMARK_DMR(
        	BISMARK_SUB.out.ch_methylation_chg,
            	BISMARK_SUB.out.ch_methylation_cpg,
            	BISMARK_SUB.out.ch_methylation_chh    
	)
        BISMARK_DMR.out.chg_bedgraph.view()
       }
    if(start =="D"){
        chg_bedgraph = tuple([[id:"group"],params.chg_bedgraph])
        cpg_bedgraph = tuple([[id:"group"],params.cpg_bedgraph])
        chh_bedgraph = tuple([[id:"group"],params.chh_bedgraph])        
        bedgraph = tuple([[id:"group"],params.chg_bedgraph,params.cpg_bedgraph,params.chh_bedgraph])       
        OSCAEWAS(
        bedgraph,
        chg_bedgraph,
        cpg_bedgraph,
        chh_bedgraph,
        PREPARE_GENOME.out.fasta_index 
     )
    }else if('D' in processList){
        OSCAEWAS(
        BISMARK_DMR.out.bedgraph,
	BISMARK_DMR.out.chg_bedgraph,
        BISMARK_DMR.out.cpg_bedgraph,
        BISMARK_DMR.out.chh_bedgraph,
        PREPARE_GENOME.out.fasta_index
    )
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
