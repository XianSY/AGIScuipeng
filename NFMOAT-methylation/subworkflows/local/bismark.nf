/*
 * bismark subworkflow
 */
include { BISMARK_ALIGN                               } from '../../modules/nf-core/bismark/align/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_ALIGNED      } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_DEDUPLICATED } from '../../modules/nf-core/samtools/sort/main'
include { BISMARK_DEDUPLICATE                         } from '../../modules/nf-core/bismark/deduplicate/main'
include { BISMARK_METHYLATIONEXTRACTOR                } from '../../modules/nf-core/bismark/methylationextractor/main'
include { BISMARK_COVERAGE2CYTOSINE                   } from '../../modules/nf-core/bismark/coverage2cytosine/main'
include { BISMARK_REPORT                              } from '../../modules/nf-core/bismark/report/main'
include { BISMARK_SUMMARY                             } from '../../modules/nf-core/bismark/summary/main'
include { BISMARK2BEDGRAPH as BISMARK2BEDGRAPH_CHG    } from '../../modules/local/bismark2graph'
include { BISMARK2BEDGRAPH as BISMARK2BEDGRAPH_CpG    } from '../../modules/local/bismark2graph'
include { BISMARK2BEDGRAPH as BISMARK2BEDGRAPH_CHH    } from '../../modules/local/bismark2graph'
include { METILENE_INPUT   as METILENE_INPUT_CpG      } from '../../modules/local/metilene_input'
include { METILENE_INPUT   as METILENE_INPUT_CHG      } from '../../modules/local/metilene_input'
include { METILENE_INPUT   as METILENE_INPUT_CHH      } from '../../modules/local/metilene_input'
include { METILENEDMR                                 } from '../../modules/local/metilenedmr'


workflow BISMARK {
    take:
    reads              // channel: [ val(meta), [ reads ] ]
    bismark_index      // channel: /path/to/BismarkIndex/
    skip_deduplication // boolean: whether to deduplicate alignments
    cytosine_report    // boolean: whether the run coverage2cytosine

    main:
    versions = Channel.empty()


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

    /*
     * Run bismark_methylation_extractor
     */
    BISMARK_METHYLATIONEXTRACTOR (
        alignments,
        bismark_index
    )
    versions = versions.mix(BISMARK_METHYLATIONEXTRACTOR.out.versions)
    //BISMARK_METHYLATIONEXTRACTOR.out.methylation_calls.view()
    
    ch_bismark2bedgraph_CHG = BISMARK_METHYLATIONEXTRACTOR.out.methylation_chg.map{
        it -> tuple(it[0],it[1],it[2])
    }.view()
    ch_bismark2bedgraph_CpG = BISMARK_METHYLATIONEXTRACTOR.out.methylation_cpg.map{
        it -> tuple(it[0],it[1],it[2])
    }
    ch_bismark2bedgraph_CHH = BISMARK_METHYLATIONEXTRACTOR.out.methylation_chh.map{
        it -> tuple(it[0],it[1],it[2])
    }
    /*
     * 分开获取得到的CpG、CHH、CHG的 beggraph 文件
     * 1、首先将CpG、CHH、CHG的calling 文件分开成三个输入文件
    */
    ch_group = channel.fromPath(params.group).map{it -> tuple(['id':'group'],it)}
    BISMARK2BEDGRAPH_CHG(
        ch_bismark2bedgraph_CHG
    )   
    BISMARK2BEDGRAPH_CpG(
        ch_bismark2bedgraph_CpG
    )
    BISMARK2BEDGRAPH_CHH(
        ch_bismark2bedgraph_CHH
    )
    //BISMARK2BEDGRAPH_CHG.out.bedgraph.map{it -> it[1]}.collect().view()
    
    METILENE_INPUT_CHG(   
    ch_group,
    BISMARK2BEDGRAPH_CHG.out.bedgraph.map{it -> it[1]}.collect()
    )
    METILENE_INPUT_CpG(
    ch_group,
    BISMARK2BEDGRAPH_CpG.out.bedgraph.map{it -> it[1]}.collect()
    )  
    METILENE_INPUT_CHH(
    ch_group,
    BISMARK2BEDGRAPH_CHH.out.bedgraph.map{it -> it[1]}.collect()
    )
    ch_compare = channel.fromPath(params.compare).splitCsv(header:true).map{row -> tuple(row.formerGroup,row.latergroup)}
    ch_dmrbedgraph = METILENE_INPUT_CHG.out.group_bedgraph.join(METILENE_INPUT_CpG.out.group_bedgraph).join(METILENE_INPUT_CHH.out.group_bedgraph).groupTuple()
    ch_dmrbedgraph.view()
    METILENEDMR(
     ch_compare,
     ch_dmrbedgraph      
    )

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
    
    emit:
    bam        = SAMTOOLS_SORT_ALIGNED.out.bam        // channel: [ val(meta), [ bam ] ] ## sorted, non-deduplicated (raw) BAM from aligner
    dedup      = SAMTOOLS_SORT_DEDUPLICATED.out.bam   // channel: [ val(meta), [ bam ] ] ## sorted, possibly deduplicated BAM
    mqc        = multiqc_files                        // path: *{html,txt}
    versions                                       // path: *.version.txt
}
