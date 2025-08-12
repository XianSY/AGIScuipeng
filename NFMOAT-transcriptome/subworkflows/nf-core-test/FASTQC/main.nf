include { INPUT_CHECK                      } from '../../../subworkflows/local/input_check'
include { CAT_FASTQ                        } from '../../../modules/nf-core/cat/fastq/main'
include { BBMAP_BBSPLIT                    } from '../../../modules/nf-core/bbmap/bbsplit/main'
include { SORTMERNA                        } from '../../../modules/nf-core/sortmerna/main'
include { FASTQC                           } from '../../../modules/nf-core/fastqc/main'
include { TRIMGALORE                       } from '../../../modules/nf-core/trimgalore/main'
include { FASTQC as FASTQC_RAW             } from '../../../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIM            } from '../../../modules/nf-core/fastqc/main'
include { UMITOOLS_EXTRACT                 } from '../../../modules/nf-core/umitools/extract/main'
include { FASTP                            } from '../../../modules/nf-core/fastp/main'
include { MULTIQC as MULTIQC_AFTER                    } from '../../../modules/nf-core/multiqc/main'
include { MULTIQC as MULTIQC_BEFORE                   } from '../../../modules/nf-core/multiqc/main'

//
// Function that parses TrimGalore log output file to get total number of reads after trimming
//
def getTrimGaloreReadsAfterFiltering(log_file) {
    def total_reads = 0
    def filtered_reads = 0
    log_file.eachLine { line ->
        def total_reads_matcher = line =~ /([\d\.]+)\ssequences processed in total/
        def filtered_reads_matcher = line =~ /shorter than the length cutoff[^:]+:\s([\d\.]+)/
        if (total_reads_matcher) total_reads = total_reads_matcher[0][1].toFloat()
        if (filtered_reads_matcher) filtered_reads = filtered_reads_matcher[0][1].toFloat()
    }
    return total_reads - filtered_reads
}



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)



//
// Function that parses fastp json output file to get total number of reads after trimming
//
import groovy.json.JsonSlurper

def getFastpReadsAfterFiltering(json_file) {
    def Map json = (Map) new JsonSlurper().parseText(json_file.text).get('summary')
    return json['after_filtering']['total_reads'].toInteger()
}

workflow FASTQC_WORKFLOW{
    take:
    ch_input
    ch_fasta
    gtf
    bbsplit_index
    prepareToolIndices
    pass_trimmed_reads

    

    main:

        ch_versions = Channel.empty()
         //
        // SUBWORKFLOW: Read in samplesheet, validate and stage input files
        //
        INPUT_CHECK (
            ch_input
        )
        .reads
        .map {
            meta, fastq ->
                new_id = meta.id - ~/_T\d+/
                [ meta + [id: new_id], fastq ]
        }
        .set { ch_fastq }
        ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

        //
        // MODULE: Concatenate FastQ files from same sample if required
        //
        CAT_FASTQ (
            ch_fastq
        )
        .reads
        .set { ch_cat_fastq }
        ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first().ifEmpty(null))

        // Branch FastQ channels if 'auto' specified to infer strandedness
        ch_cat_fastq
            .branch {
                meta, fastq ->
                    auto_strand : meta.strandedness == 'auto'
                        return [ meta, fastq ]
                    known_strand: meta.strandedness != 'auto'
                        return [ meta, fastq ]
            }
            .set { ch_strand_fastq }
          
        fastqc_html = Channel.empty()
        fastqc_zip  = Channel.empty()
	FASTQC_RAW(ch_strand_fastq.auto_strand)
        fastqc_html = FASTQC_RAW.out.html
        fastqc_zip  = FASTQC_RAW.out.zip
	ch_versions = ch_versions.mix(FASTQC_RAW.out.versions.first())
 
        ch_multiqc_files = Channel.empty()
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC_RAW.out.zip.collect{it[1]}.ifEmpty([]))
        MULTIQC_BEFORE (
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList()
        )


        //
        // SUBWORKFLOW: Read QC, extract UMI and trim adapters with TrimGalore!
        //
        ch_filtered_reads      = Channel.empty()
        ch_fastqc_raw_multiqc  = Channel.empty()
        ch_fastqc_trim_multiqc = Channel.empty()
        ch_trim_log_multiqc    = Channel.empty()
        ch_trim_read_count     = Channel.empty()
        if (params.trimmer == 'trimgalore') {

            umi_reads = ch_strand_fastq.auto_strand
            umi_log   = Channel.empty()
            if (params.with_umi && !params.skip_umi_extract) {
                UMITOOLS_EXTRACT (ch_strand_fastq)
                umi_reads   = UMITOOLS_EXTRACT.out.reads
                umi_log     = UMITOOLS_EXTRACT.out.log
                ch_versions = ch_versions.mix(UMITOOLS_EXTRACT.out.versions.first())

                // Discard R1 / R2 if required
                if (params.umi_discard_read in [1,2]) {
                    UMITOOLS_EXTRACT
                        .out
                        .reads
                        .map {
                            meta, reads ->
                                meta.single_end ? [ meta, reads ] : [ meta + ['single_end': true], reads[params.umi_discard_read % 2] ]
                        }
                        .set { umi_reads }
                }
            }
        trim_reads      = umi_reads
        trim_unpaired   = Channel.empty()
        trim_html       = Channel.empty()
        trim_zip        = Channel.empty()
        trim_log        = Channel.empty()
        trim_read_count = Channel.empty()
        if (!params.skip_trimming) {
            TRIMGALORE (umi_reads)
            trim_unpaired = TRIMGALORE.out.unpaired
            trim_html     = TRIMGALORE.out.html
            trim_zip      = TRIMGALORE.out.zip
            trim_log      = TRIMGALORE.out.log
            ch_versions   = ch_versions.mix(TRIMGALORE.out.versions.first())

            //
            // Filter FastQ files based on minimum trimmed read count after adapter trimming
            //
            TRIMGALORE
                .out
                .reads
                .join(trim_log, remainder: true)
                .map {
                    meta, reads, trim_log ->
                        if (trim_log) {
                            num_reads = getTrimGaloreReadsAfterFiltering(meta.single_end ? trim_log : trim_log[-1])
                            [ meta, reads, num_reads ]
                        } else {
                            [ meta, reads, params.min_trimmed_reads.toFloat() + 1 ]
                        }
                }
                .set { ch_num_trimmed_reads }

                ch_num_trimmed_reads
                    .filter { meta, reads, num_reads -> num_reads >= params.min_trimmed_reads.toFloat() }
                    .map { meta, reads, num_reads -> [ meta, reads ] }
                    .set { trim_reads }

                ch_num_trimmed_reads
                    .map { meta, reads, num_reads -> [ meta, num_reads ] }
                    .set { trim_read_count }
            }
            
            FASTQC_TRIM(trim_reads)
            ch_filtered_reads      = trim_reads
            ch_fastqc_raw_multiqc  = fastqc_html
            ch_fastqc_trim_multiqc = trim_zip
            ch_trim_log_multiqc    = umi_log
            ch_trim_read_count     = trim_read_count
	    ch_multiqc_zip = FASTQC_TRIM.out.zip.collect{it[1]}.ifEmpty([])
        }

        //
        // SUBWORKFLOW: Read QC, extract UMI and trim adapters with fastp
        //
        if (params.trimmer == 'fastp') {

            umi_reads = reads
            umi_log   = Channel.empty()
            if (params.with_umi && !params.skip_umi_extract) {
                UMITOOLS_EXTRACT (
                    ch_strand_fastq.auto_strand
                )
                umi_reads   = UMITOOLS_EXTRACT.out.reads
                umi_log     = UMITOOLS_EXTRACT.out.log
                ch_versions = ch_versions.mix(UMITOOLS_EXTRACT.out.versions.first())

                // Discard R1 / R2 if required
                if (params.umi_discard_read in [1,2]) {
                    UMITOOLS_EXTRACT
                        .out
                        .reads
                        .map {
                            meta, reads ->
                                meta.single_end ? [ meta, reads ] : [ meta + [single_end: true], reads[params.umi_discard_read % 2] ]
                        }
                        .set { umi_reads }
                }
            }

            trim_reads        = umi_reads
            trim_json         = Channel.empty()
            trim_html         = Channel.empty()
            trim_log          = Channel.empty()
            trim_reads_fail   = Channel.empty()
            trim_reads_merged = Channel.empty()
            fastqc_trim_html  = Channel.empty()
            fastqc_trim_zip   = Channel.empty()
            trim_read_count   = Channel.empty()

            if (!params.skip_trimming) {
                FASTP (
                    umi_reads,
                    [],
                    params.save_trimmed,
                    params.save_trimmed
                )
                trim_json         = FASTP.out.json
                trim_html         = FASTP.out.html
                trim_log          = FASTP.out.log
                trim_reads_fail   = FASTP.out.reads_fail
                trim_reads_merged = FASTP.out.reads_merged
                ch_versions       = ch_versions.mix(FASTP.out.versions.first())

            //
            // Filter FastQ files based on minimum trimmed read count after adapter trimming
            //
            FASTP
                .out
                .reads
                .join(trim_json)
                .map { meta, reads, json -> [ meta, reads, getFastpReadsAfterFiltering(json) ] }
                .set { ch_num_trimmed_reads }

            ch_num_trimmed_reads
                .filter { meta, reads, num_reads -> num_reads >= params.min_trimmed_reads.toInteger() }
                .map { meta, reads, num_reads -> [ meta, reads ] }
                .set { trim_reads }

            ch_num_trimmed_reads
                .map { meta, reads, num_reads -> [ meta, num_reads ] }
                .set { trim_read_count }

            FASTQC_TRIM (
                    trim_reads
                )
                fastqc_trim_html = FASTQC_TRIM.out.html
                fastqc_trim_zip  = FASTQC_TRIM.out.zip
                ch_versions      = ch_versions.mix(FASTQC_TRIM.out.versions.first())
            
        }

            ch_filtered_reads      = trim_reads
            ch_fastqc_raw_multiqc  = fastqc_raw_zip
            ch_fastqc_trim_multiqc = fastqc_trim_zip
            ch_trim_log_multiqc    = trim_json
            ch_trim_read_count     = trim_read_count
            ch_multiqc_zip = FASTQC_TRIM.out.zip.collect{it[1]}.ifEmpty([])
       }

         
        ch_multiqc_files = Channel.empty()
        ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_zip)

        MULTIQC_AFTER (
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList()
        )
 
        //
        // Get list of samples that failed trimming threshold for MultiQC report
        //
        ch_trim_read_count
            .map {
                meta, num_reads ->
                    pass_trimmed_reads[meta.id] = true
                    if (num_reads <= params.min_trimmed_reads.toFloat()) {
                        pass_trimmed_reads[meta.id] = false
                        return [ "$meta.id\t$num_reads" ]
                    }
            }
            .collect()
            .map {
                tsv_data ->
                    def header = ["Sample", "Reads after trimming"]
                    WorkflowRnaseq.multiqcTsvFromList(tsv_data, header)
            }
            .set { ch_fail_trimming_multiqc }

        //
        // MODULE: Remove genome contaminant reads
        //
        if (!params.skip_bbsplit) {
            BBMAP_BBSPLIT (
                ch_filtered_reads,
                PREPARE_GENOME.out.bbsplit_index,
                [],
                [ [], [] ],
                false
            )
            .primary_fastq
            .set { ch_filtered_reads }
            ch_versions = ch_versions.mix(BBMAP_BBSPLIT.out.versions.first())
        }

        //
        // MODULE: Remove ribosomal RNA reads
        //
        ch_sortmerna_multiqc = Channel.empty()
        if (params.remove_ribo_rna) {
            ch_sortmerna_fastas = Channel.from(ch_ribo_db.readLines()).map { row -> file(row, checkIfExists: true) }.collect()

            SORTMERNA (
                ch_filtered_reads,
                ch_sortmerna_fastas
            )
            .reads
            .set { ch_filtered_reads }

            ch_sortmerna_multiqc = SORTMERNA.out.log
            ch_versions = ch_versions.mix(SORTMERNA.out.versions.first())
        }
    emit:
    ch_filtered_reads
    ch_versions 
       
}
