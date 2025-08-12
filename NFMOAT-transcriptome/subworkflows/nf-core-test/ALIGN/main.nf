include { ALIGN_STAR                       } from '../../../subworkflows/local/align_star'
include { SAMTOOLS_SORT                    } from '../../../modules/nf-core/samtools/sort/main'
include { UMITOOLS_PREPAREFORRSEM as UMITOOLS_PREPAREFORSALMON } from '../../../modules/local/umitools_prepareforrsem.nf'
include { UMITOOLS_DEDUP     as UMITOOLS_DEDUP_GENOME            } from '../../../modules/nf-core/umitools/dedup/main'
include { SAMTOOLS_INDEX     as SAMTOOLS_INDEX_GENOME            } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_STATS    as SAMTOOLS_STATS_GENOME             } from '../../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_IDXSTATS as SAMTOOLS_IDXSTATS_GENOME          } from '../../../modules/nf-core/samtools/idxstats/main'
include { SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_GENOME          } from '../../../modules/nf-core/samtools/flagstat/main'

include { UMITOOLS_DEDUP     as UMITOOLS_DEDUP_TRANSCRIPTOME     } from '../../../modules/nf-core/umitools/dedup/main'
include { SAMTOOLS_INDEX     as SAMTOOLS_INDEX_TRANSCRIPTOME     } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_STATS    as SAMTOOLS_STATS_TRANSCRIPTOME      } from '../../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_IDXSTATS as SAMTOOLS_IDXSTATS_TRANSCRIPTOME   } from '../../../modules/nf-core/samtools/idxstats/main'
include { SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_TRANSCRIPTOME   } from '../../../modules/nf-core/samtools/flagstat/main'
//include { SAMTOOLS_SORT                                          } from '../../../modules/nf-core/samtools/sort/main'

include { SAMTOOLS_SORT      as SAMTOOLS_SORT_TRANSCRIPTOME      } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX     } from '../../../modules/nf-core/samtools/index/main'

workflow ALIGN{
        take:
        ch_filtered_reads
        star_index
        gtf
        fasta_map
        hisat2_index_map
        splicesites
        ch_versions
        is_aws_igenome

        main:
            //
            // SUBWORKFLOW: Alignment with STAR and gene/transcript quantification with Salmon
            //
            ch_genome_bam                 = Channel.empty()
            ch_genome_bam_index           = Channel.empty()
            ch_samtools_stats             = Channel.empty()
            ch_samtools_flagstat          = Channel.empty()
            ch_samtools_idxstats          = Channel.empty()
            //
            ch_genome_bam                 = Channel.empty()
            ch_genome_bam_index           = Channel.empty()
            ch_samtools_stats             = Channel.empty()
            ch_samtools_flagstat          = Channel.empty()
            ch_samtools_idxstats          = Channel.empty()
            ch_star_multiqc               = Channel.empty()
            ch_aligner_pca_multiqc        = Channel.empty()
            ch_aligner_clustering_multiqc = Channel.empty()

            // align with star tools
            if (!params.skip_alignment && params.aligner == 'star_salmon') {
                ALIGN_STAR (
                    ch_filtered_reads,
                    star_index,
                    gtf,
                    params.star_ignore_sjdbgtf,
                    '',
                    params.seq_center ?: '',
                    is_aws_igenome,
                    fasta_map
                )
            ch_genome_bam        = ALIGN_STAR.out.bam
            ch_genome_bam_index  = ALIGN_STAR.out.bai
            ch_transcriptome_bam = ALIGN_STAR.out.bam_transcript
            ch_samtools_stats    = ALIGN_STAR.out.stats
            ch_samtools_flagstat = ALIGN_STAR.out.flagstat
            ch_samtools_idxstats = ALIGN_STAR.out.idxstats
            ch_star_multiqc      = ALIGN_STAR.out.log_final
            ch_genome_bam_total  = ALIGN_STAR.out.ch_genome_bam_total
            if (params.bam_csi_index) {
                ch_genome_bam_index = ALIGN_STAR.out.csi
                }
            ch_versions = ch_versions.mix(ALIGN_STAR.out.versions)
            }

            //
        // SUBWORKFLOW: Remove duplicate reads from BAM file based on UMIs
        // 去重
        if (params.with_umi) { //default with_umi = false
            // Deduplicate genome BAM file before downstream analysis

            // from bam_dedup_stats_samtools_umitools_GENOME
            // UMI-tools dedup
            //

        

        UMITOOLS_DEDUP_GENOME (
                ch_genome_bam.join(ch_genome_bam_index, by: [0]),
                params.umitools_dedup_stats
             )
            ch_versions = ch_versions.mix(UMITOOLS_DEDUP_GENOME.out.versions.first())

            //
            // Index BAM file and run samtools stats, flagstat and idxstats
            //

            SAMTOOLS_INDEX_GENOME ( UMITOOLS_DEDUP_GENOME.out.bam )
            ch_versions = ch_versions.mix(SAMTOOLS_INDEX_GENOME.out.versions.first())

            ch_bam_bai_dedup = UMITOOLS_DEDUP_GENOME.out.bam
                .join(SAMTOOLS_INDEX_GENOME.out.bai, by: [0], remainder: true)
                .join(SAMTOOLS_INDEX_GENOME.out.csi, by: [0], remainder: true)
                .map {
                    meta, bam, bai, csi ->
                        if (bai) {
                            [ meta, bam, bai ]
                        } else {
                            [ meta, bam, csi ]
                        }
                }

            SAMTOOLS_STATS_GENOME ( ch_bam_bai_dedup, [ [:], [] ] )
            ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions)

            SAMTOOLS_FLAGSTAT_GENOME ( ch_bam_bai_dedup )
            ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions)

            SAMTOOLS_IDXSTATS_GENOME ( ch_bam_bai_dedup )
            ch_versions = ch_versions.mix(SAMTOOLS_IDXSTATS.out.versions)

            ch_genome_bam        = UMITOOLS_DEDUP_GENOME.out.bam
            ch_genome_bam_index  = SAMTOOLS_INDEX_GENOME.out.bai
            ch_samtools_stats    = SAMTOOLS_STATS_GENOME.out.stats
            ch_samtools_flagstat = SAMTOOLS_FLAGSTAT_GENOME.out.flagstat
            ch_samtools_idxstats = SAMTOOLS_IDXSTATS_GENOME.out.idxstats

            if (params.bam_csi_index) {
                ch_genome_bam_index  = SAMTOOLS_INDEX_GENOME.out.csi
            }
            // Co-ordinate sort, index and run stats on transcriptome BAM
            //from bam_sort_stats_samtools

            SAMTOOLS_SORT_TRANSCRIPTOME ( ch_transcriptome_bam )
            ch_versions = ch_versions.mix(SAMTOOLS_SORT_TRANSCRIPTOME.out.versions.first())

            SAMTOOLS_INDEX_TRANSCRIPTOME ( SAMTOOLS_SORT_TRANSCRIPTOME.out.bam )
            ch_versions = ch_versions.mix(SAMTOOLS_INDEX_TRANSCRIPTOME.out.versions.first())

            SAMTOOLS_SORT_TRANSCRIPTOME.out.bam
                .join(SAMTOOLS_INDEX_TRANSCRIPTOME.out.bai, by: [0], remainder: true)
                .join(SAMTOOLS_INDEX_TRANSCRIPTOME.out.csi, by: [0], remainder: true)
                .map {
                    meta, bam, bai, csi ->
                        if (bai) {
                            [ meta, bam, bai ]
                        } else {
                            [ meta, bam, csi ]
                        }
                }
                .set { ch_bam_bai }        

            BAM_STATS_SAMTOOLS_TRANSCRIPTOME ( ch_bam_bai, fasta_map )
            ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS_TRANSCRIPTOME.out.versions)

            ch_transcriptome_sorted_bam = SAMTOOLS_SORT_TRANSCRIPTOME.out.bam
            ch_transcriptome_sorted_bai = SAMTOOLS_INDEX_TRANSCRIPTOME.out.bai


            // Deduplicate transcriptome BAM file before read counting with Salmon

            // from bam_dedup_stats_samtools_umitools_transcriptome
            // UMI-tools dedup
            //
            UMITOOLS_DEDUP_TRANSCRIPTOME (
                ch_transcriptome_sorted_bam.join(ch_transcriptome_sorted_bai, by: [0]),
                params.umitools_dedup_stats
            )
            ch_versions = ch_versions.mix(UMITOOLS_DEDUP_TRANSCRIPTOME.out.versions.first())

            //
            // Index BAM file and run samtools stats, flagstat and idxstats
            //
            SAMTOOLS_INDEX_TRANSCRIPTOME (UMITOOLS_DEDUP_TRANSCRIPTOME.out.bam )
            ch_versions = ch_versions.mix(SAMTOOLS_INDEX_TRANSCRIPTOME.out.versions.first())
            
            ch_bam_bai_dedup_transcript = UMITOOLS_DEDUP_TRANSCRIPTOME.out.bam
                .join(SAMTOOLS_INDEX_TRANSCRIPTOME.out.bai, by: [0], remainder: true)
                .join(SAMTOOLS_INDEX_TRANSCRIPTOME.out.csi, by: [0], remainder: true)
                .map {
                    meta, bam, bai, csi ->
                        if (bai) {
                            [ meta, bam, bai ]
                        } else {
                            [ meta, bam, csi ]
                        }
                }

            SAMTOOLS_STATS_TRANSCRIPTOME ( ch_bam_bai_dedup_transcript, [ [:], [] ] )
            ch_versions = ch_versions.mix(SAMTOOLS_STATS_TRANSCRIPTOME.out.versions)

            SAMTOOLS_FLAGSTAT_TRANSCRIPTOME ( ch_bam_bai_dedup_transcript )
            ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT_TRANSCRIPTOME.out.versions)

            SAMTOOLS_IDXSTATS_TRANSCRIPTOME ( ch_bam_bai_dedup_transcript )
            ch_versions = ch_versions.mix(SAMTOOLS_IDXSTATS_TRANSCRIPTOME.out.versions)


            // Name sort BAM before passing to Salmon
            SAMTOOLS_SORT (
                UMITOOLS_DEDUP_TRANSCRIPTOME.out.bam
            )

            
            // Name sort BAM before passing to Salmon
            SAMTOOLS_SORT (
                UMITOOLS_DEDUP_TRANSCRIPTOME.out.bam
            )

            // Only run prepare_for_rsem.py on paired-end BAM files
            SAMTOOLS_SORT
                .out
                .bam
                .branch {
                    meta, bam ->
                        single_end: meta.single_end
                            return [ meta, bam ]
                        paired_end: !meta.single_end
                            return [ meta, bam ]
                }
                .set { ch_umitools_dedup_bam_transcript }

            // Fix paired-end reads in name sorted BAM file
            // See: https://github.com/nf-core/rnaseq/issues/828
            UMITOOLS_PREPAREFORSALMON (
                ch_umitools_dedup_bam_transcript.paired_end
            )
            ch_versions = ch_versions.mix(UMITOOLS_PREPAREFORSALMON.out.versions.first())

            ch_umitools_dedup_bam
                .single_end
                .mix(UMITOOLS_PREPAREFORSALMON.out.bam)
                .set { ch_transcriptome_bam }
        }

        // hisat2比对

        //
        // SUBWORKFLOW: Alignment with HISAT2
        //
        ch_hisat2_multiqc = Channel.empty()
        if (!params.skip_alignment && params.aligner == 'hisat2') {
            FASTQ_ALIGN_HISAT2 (
                ch_filtered_reads,
                PREPARE_GENOME.out.hisat2_index.map { [ [:], it ] },
                PREPARE_GENOME.out.splicesites.map { [ [:], it ] },
                PREPARE_GENOME.out.fasta.map { [ [:], it ] }
            )
            ch_genome_bam        = FASTQ_ALIGN_HISAT2.out.bam
            ch_genome_bam_index  = FASTQ_ALIGN_HISAT2.out.bai
            ch_samtools_stats    = FASTQ_ALIGN_HISAT2.out.stats
            ch_samtools_flagstat = FASTQ_ALIGN_HISAT2.out.flagstat
            ch_samtools_idxstats = FASTQ_ALIGN_HISAT2.out.idxstats
            ch_hisat2_multiqc    = FASTQ_ALIGN_HISAT2.out.summary
            if (params.bam_csi_index) {
                ch_genome_bam_index = FASTQ_ALIGN_HISAT2.out.csi
            }
            ch_versions = ch_versions.mix(FASTQ_ALIGN_HISAT2.out.versions)
            //
            // SUBWORKFLOW: Remove duplicate reads from BAM file based on UMIs
            //
            if (params.with_umi) {
                BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME (
                    ch_genome_bam.join(ch_genome_bam_index, by: [0]),
                    params.umitools_dedup_stats
                )
                ch_genome_bam        = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.bam
                ch_genome_bam_index  = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.bai
                ch_samtools_stats    = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.stats
                ch_samtools_flagstat = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.flagstat
                ch_samtools_idxstats = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.idxstats
                ch_genome_bam_total = BAM_MARKDUPLICATES_PICARD.out.bam.map{it->it[1]}.collect()
                if (params.bam_csi_index) {
                    ch_genome_bam_index = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.csi
                }
                ch_versions = ch_versions.mix(BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.versions)
            }


        }      

    emit:
        ch_genome_bam
        ch_genome_bam_total
}


