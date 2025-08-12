//
// Alignment with STAR
//

include { STAR_ALIGN          } from '../../modules/nf-core/star/align/main'
include { STAR_ALIGN_IGENOMES } from '../../modules/local/star_align_igenomes'
include { SAMTOOLS_SORT      } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX     } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_STATS    } from '../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_IDXSTATS } from '../../modules/nf-core/samtools/idxstats/main'
include { SAMTOOLS_FLAGSTAT } from '../../modules/nf-core/samtools/flagstat/main'


workflow ALIGN_STAR {
    take:
    reads               // channel: [ val(meta), [ reads ] ]
    index               // channel: /path/to/star/index/
    gtf                 // channel: /path/to/genome.gtf
    star_ignore_sjdbgtf // boolean: when using pre-built STAR indices do not re-extract and use splice junctions from the GTF file 
    seq_platform        // string : sequencing platform
    seq_center          // string : sequencing center
    is_aws_igenome      // boolean: whether the genome files are from AWS iGenomes
    fasta               // channel: /path/to/fasta

    main:

    ch_versions = Channel.empty()

    //
    // Map reads with STAR
    //
    ch_orig_bam       = Channel.empty()
    ch_log_final      = Channel.empty()
    ch_log_out        = Channel.empty()
    ch_log_progress   = Channel.empty()
    ch_bam_sorted     = Channel.empty()
    ch_bam_transcript = Channel.empty()
    ch_fastq          = Channel.empty()
    ch_tab            = Channel.empty()
    if (is_aws_igenome) {
        STAR_ALIGN_IGENOMES ( reads, index, gtf, star_ignore_sjdbgtf, seq_platform, seq_center )
        ch_orig_bam       = STAR_ALIGN_IGENOMES.out.bam
        ch_log_final      = STAR_ALIGN_IGENOMES.out.log_final
        ch_log_out        = STAR_ALIGN_IGENOMES.out.log_out
        ch_log_progress   = STAR_ALIGN_IGENOMES.out.log_progress
        ch_bam_sorted     = STAR_ALIGN_IGENOMES.out.bam_sorted
        ch_bam_transcript = STAR_ALIGN_IGENOMES.out.bam_transcript
        ch_fastq          = STAR_ALIGN_IGENOMES.out.fastq
        ch_tab            = STAR_ALIGN_IGENOMES.out.tab
        ch_versions       = ch_versions.mix(STAR_ALIGN_IGENOMES.out.versions.first())
    } else {
        STAR_ALIGN ( reads, index, gtf, star_ignore_sjdbgtf, seq_platform, seq_center )
        ch_orig_bam       = STAR_ALIGN.out.bam
        ch_log_final      = STAR_ALIGN.out.log_final
        ch_log_out        = STAR_ALIGN.out.log_out
        ch_log_progress   = STAR_ALIGN.out.log_progress
        ch_bam_sorted     = STAR_ALIGN.out.bam_sorted
        ch_bam_transcript = STAR_ALIGN.out.bam_transcript
        ch_fastq          = STAR_ALIGN.out.fastq
        ch_tab            = STAR_ALIGN.out.tab
        ch_versions       = ch_versions.mix(STAR_ALIGN.out.versions.first())
    }

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    SAMTOOLS_SORT ( ch_orig_bam )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    SAMTOOLS_INDEX ( SAMTOOLS_SORT.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    SAMTOOLS_SORT.out.bam
        .join(SAMTOOLS_INDEX.out.bai, by: [0], remainder: true)
        .join(SAMTOOLS_INDEX.out.csi, by: [0], remainder: true)
        .map {
            meta, bam, bai, csi ->
                if (bai) {
                    [ meta, bam, bai ]
                } else {
                    [ meta, bam, csi ]
                }
        }
        .set { ch_bam_bai }

    SAMTOOLS_STATS ( ch_bam_bai, fasta )
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions)

    SAMTOOLS_FLAGSTAT ( ch_bam_bai )
    ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions)

    SAMTOOLS_IDXSTATS ( ch_bam_bai )
    ch_versions = ch_versions.mix(SAMTOOLS_IDXSTATS.out.versions)


    ch_genome_bam_total = SAMTOOLS_SORT.out.bam.map{it->it[1]}.collect()

    emit:
    orig_bam       = ch_orig_bam                    // channel: [ val(meta), bam            ]
    log_final      = ch_log_final                   // channel: [ val(meta), log_final      ]
    log_out        = ch_log_out                     // channel: [ val(meta), log_out        ]
    log_progress   = ch_log_progress                // channel: [ val(meta), log_progress   ]
    bam_sorted     = ch_bam_sorted                  // channel: [ val(meta), bam_sorted     ]
    bam_transcript = ch_bam_transcript              // channel: [ val(meta), bam_transcript ]
    fastq          = ch_fastq                       // channel: [ val(meta), fastq          ]
    tab            = ch_tab                         // channel: [ val(meta), tab            ]

    bam            = SAMTOOLS_SORT.out.bam      // channel: [ val(meta), [ bam ] ]
    bai            = SAMTOOLS_INDEX.out.bai      // channel: [ val(meta), [ bai ] ]
    csi            = SAMTOOLS_INDEX.out.csi      // channel: [ val(meta), [ csi ] ]
    stats          = SAMTOOLS_STATS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat       = SAMTOOLS_FLAGSTAT.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats       = SAMTOOLS_IDXSTATS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
    ch_genome_bam_total                             //channel: [ val(meta),[bam1,bam2,..]]

    versions       = ch_versions                    // channel: [ versions.yml ]
}
