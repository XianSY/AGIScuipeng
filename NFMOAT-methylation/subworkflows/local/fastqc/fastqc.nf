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



include { CAT_FASTQ                   } from '../../../modules/nf-core/cat/fastq/main'
include { FASTQC  as FASTQC_BEFORE    } from '../../../modules/nf-core/fastqc/main'
include { FASTQC  as FASTQC_AFTER     } from '../../../modules/nf-core/fastqc/main'
include { MULTIQC as MULTIQC_BEFRORE  } from '../../../modules/nf-core/multiqc/main'
include { MULTIQC as MULTIQC_CLEAN    } from '../../../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../../../modules/nf-core/custom/dumpsoftwareversions/main'
include { TRIMGALORE                  } from '../../../modules/nf-core/trimgalore/main'


workflow FASTQC_WORK{
    main:
       ch_versions = Channel.empty()
       channel.fromPath(params.input)
        .splitCsv(header:true)
        .map {
            data ->
            if (! data.fastq_2 ) {
                return [ ["id":data.sample, "single_end":true], [ data.fastq_1 ] ]
            } else {
                return [ ["id":data.sample, "single_end":false], [data.fastq_1, data.fastq_2 ] ]
            }
        }.set { ch_fastq }
        //
        // MODULE: Concatenate FastQ files from same sample if required
        //
        CAT_FASTQ (
            ch_fastq
        )
        .reads
        .set { ch_cat_fastq }
        ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first())
        

        //
        // MODULE: Run FastQC
        //
        FASTQC_BEFORE (
            ch_cat_fastq
        )
        ch_versions = ch_versions.mix(FASTQC_BEFORE.out.versions.first())

        
        CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
        )

        //
        // MODULE: MultiQC
        //
        //WorkflowTest.groovy
        workflow_summary    = WorkflowMethylseq.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        methods_description    = WorkflowMethylseq.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
        ch_methods_description = Channel.value(methods_description)

        ch_multiqc_files = Channel.empty()
        ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC_BEFORE.out.zip.collect{it[1]}.ifEmpty([]))

        MULTIQC_BEFRORE (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
        )



        /*
        * MODULE: Run TrimGalore!
        */
        if (!params.skip_trimming) {
            TRIMGALORE(ch_cat_fastq)
            reads = TRIMGALORE.out.reads
            ch_versions = ch_versions.mix(TRIMGALORE.out.versions.first())
        } else {
            reads = ch_cat_fastq
        }
        FASTQC_AFTER(reads)
        ch_versions = ch_versions.mix(FASTQC_AFTER.out.versions.first())
        ch_multiqc_files = Channel.empty()
            ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
            ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
            ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
            ch_multiqc_files = ch_multiqc_files.mix(FASTQC_AFTER.out.zip.collect{it[1]}.ifEmpty([]))
            MULTIQC_CLEAN(
                ch_multiqc_files.collect(),
                ch_multiqc_config.toList(),
                ch_multiqc_custom_config.toList(),
                ch_multiqc_logo.toList()
            )
 

     emit:
      reads = reads
}
