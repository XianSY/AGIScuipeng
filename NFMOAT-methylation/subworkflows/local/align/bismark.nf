
// Aligner: bismark or bismark_hisat
if( params.aligner =~ /bismark/ ){
    include { BISMARK } from '../bismark'
}
// Aligner: bwameth
else if ( params.aligner == 'bwameth' ){
    include { BWAMETH } from '../bwameth'
}

/*
 * bismark modules of bismark subworkflow 
 */
include { BISMARK_ALIGN                               } from '../../../modules/nf-core/bismark/align/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_ALIGNED      } from '../../../modules/nf-core/samtools/sort/main'
include { BISMARK_DEDUPLICATE                         } from '../../../modules/nf-core/bismark/deduplicate/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_DEDUPLICATED } from '../../../modules/nf-core/samtools/sort/main'
include { BISMARK_METHYLATIONEXTRACTOR                } from '../../../modules/nf-core/bismark/methylationextractor/main'
include { BISMARK_COVERAGE2CYTOSINE                   } from '../../../modules/nf-core/bismark/coverage2cytosine/main'
include { BISMARK_REPORT                              } from '../../../modules/nf-core/bismark/report/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                 } from '../../../modules/nf-core/custom/dumpsoftwareversions/main'
include { QUALIMAP_BAMQC                              } from '../../../modules/nf-core/qualimap/bamqc/main'
include { PRESEQ_LCEXTRAP                             } from '../../../modules/nf-core/preseq/lcextrap/main'
include { BISMARK_SUMMARY                             } from '../../../modules/nf-core/bismark/summary/main'
include { MULTIQC                                     } from '../../../modules/nf-core/multiqc/main'


workflow BISMARK_SUB{
        take:
        start
        reads              // channel: [ val(meta), [ reads ] ]
        bismark_index      // channel: /path/to/BismarkIndex/
        skip_deduplication // boolean: whether to deduplicate alignments
        cytosine_report    // boolean: whether the run coverage2cytosine
        summary_params
        ch_multiqc_custom_methods_description
        ch_multiqc_config
        ch_multiqc_custom_config
        ch_multiqc_logo

        main:
        versions = Channel.empty()
        ch_versions = Channel.empty()
        
        /*
        * SUBWORKFLOW: Align reads, deduplicate and extract methylation with Bismark
        */
        // Aligner: bismark or bismark_hisat
        if( params.aligner =~ /bismark/ ){
            /*
            * Run Bismark alignment + downstream processing
            */
            /*
            * Align with bismark
            */
            BISMARK_ALIGN (
                reads,
                bismark_index
            )
            versions = versions.mix(BISMARK_ALIGN.out.versions)

            /*
            * Sort raw output BAM
            */
            SAMTOOLS_SORT_ALIGNED(
                BISMARK_ALIGN.out.bam,
            )
            versions = versions.mix(SAMTOOLS_SORT_ALIGNED.out.versions)

            if (skip_deduplication) {
                alignments = BISMARK_ALIGN.out.bam
                alignment_reports = BISMARK_ALIGN.out.report.map{ meta, report -> [ meta, report, [] ] }
            } else {
                /*
                * Run deduplicate_bismark
                */
                BISMARK_DEDUPLICATE( BISMARK_ALIGN.out.bam )

                alignments = BISMARK_DEDUPLICATE.out.bam
                alignment_reports = BISMARK_ALIGN.out.report.join(BISMARK_DEDUPLICATE.out.report)
                versions = versions.mix(BISMARK_DEDUPLICATE.out.versions)
            }
            //alignment_reports.view()
            /*
            * Run bismark_methylation_extractor
            */
            BISMARK_METHYLATIONEXTRACTOR (
                alignments,
                bismark_index
            )
            versions = versions.mix(BISMARK_METHYLATIONEXTRACTOR.out.versions)
            //BISMARK_METHYLATIONEXTRACTOR.out.methylation_calls.view()
            
            /*
            * Run coverage2cytosine
            */
            if (cytosine_report) {
                BISMARK_COVERAGE2CYTOSINE (
                    BISMARK_METHYLATIONEXTRACTOR.out.coverage,
                    bismark_index
                )
                versions = versions.mix(BISMARK_COVERAGE2CYTOSINE.out.versions)
            }

            /*
            * Generate bismark sample reports
            */
            BISMARK_REPORT (
                alignment_reports
                    .join(BISMARK_METHYLATIONEXTRACTOR.out.report)
                    .join(BISMARK_METHYLATIONEXTRACTOR.out.mbias)
            )
            versions = versions.mix(BISMARK_REPORT.out.versions)

            /*
            * Generate bismark summary report
            */
            BISMARK_SUMMARY (
              BISMARK_ALIGN.out.bam.collect{ it[1].name }.ifEmpty([]),
              alignment_reports.collect{ it[1] }.ifEmpty([]),
              alignment_reports.collect{ it[2] }.ifEmpty([]),
              BISMARK_METHYLATIONEXTRACTOR.out.report.collect{ it[1] }.ifEmpty([]),
              BISMARK_METHYLATIONEXTRACTOR.out.mbias.collect{ it[1] }.ifEmpty([])
            )
            versions = versions.mix(BISMARK_SUMMARY.out.versions)
	    
            /*
            * MODULE: Run samtools sort
            */
            SAMTOOLS_SORT_DEDUPLICATED (
                alignments
            )
            versions = versions.mix(SAMTOOLS_SORT_DEDUPLICATED.out.versions)

            /*
            * Collect MultiQC inputs
            */
            BISMARK_SUMMARY.out.summary.ifEmpty([])
               .mix(alignment_reports.collect{ it[1] })
                .mix(alignment_reports.collect{ it[2] })
                .mix(BISMARK_METHYLATIONEXTRACTOR.out.report.collect{ it[1] })
                .mix(BISMARK_METHYLATIONEXTRACTOR.out.mbias.collect{ it[1] })
                .mix(BISMARK_REPORT.out.report.collect{ it[1] })
                .set{ multiqc_files }
            
             ch_aligner_mqc = multiqc_files
    
    /*
     * MODULE: Qualimap BamQC
     */
    QUALIMAP_BAMQC (
        SAMTOOLS_SORT_DEDUPLICATED.out.bam,
        params.bamqc_regions_file ? Channel.fromPath( params.bamqc_regions_file, checkIfExists: true ).toList() : []
    )
    ch_versions = ch_versions.mix(QUALIMAP_BAMQC.out.versions.first())

    /*
     * MODULE: Run Preseq
     */
    PRESEQ_LCEXTRAP (
        SAMTOOLS_SORT_ALIGNED.out.bam
    )
    ch_versions = ch_versions.mix(PRESEQ_LCEXTRAP.out.versions.first())

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    if (!params.skip_multiqc) {
        workflow_summary    = WorkflowMethylseq.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        methods_description    = WorkflowMethylseq.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
        ch_methods_description = Channel.value(methods_description)

        ch_multiqc_files = Channel.empty()
        ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
        ch_multiqc_files = ch_multiqc_files.mix(QUALIMAP_BAMQC.out.results.collect{ it[1] }.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(PRESEQ_LCEXTRAP.out.log.collect{ it[1] }.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_aligner_mqc.ifEmpty([]))
        
        if (start == "E" & !params.skip_trimming) {
            ch_multiqc_files = ch_multiqc_files.mix(log.collect{ it[1] })
            ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{ it[1] }.ifEmpty([]))
        }
        

        MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList()
        )
        multiqc_report = MULTIQC.out.report.toList()
        ch_versions    = ch_versions.mix(MULTIQC.out.versions)
       //BISMARK_METHYLATIONEXTRACTOR.out.methylation_chg.view()
     }
        
   }
   emit:
     ch_methylation_chg = BISMARK_METHYLATIONEXTRACTOR.out.methylation_chg
     ch_methylation_cpg = BISMARK_METHYLATIONEXTRACTOR.out.methylation_cpg
     ch_methylation_chh = BISMARK_METHYLATIONEXTRACTOR.out.methylation_chh
     //versions                                       // path: *.version.txt
   
}
