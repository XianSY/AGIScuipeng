include { GATK4_CREATESEQUENCEDICTIONARY                           } from '../../../modules/nf-core/gatk4/createsequencedictionary/main'
include { GATK4_HAPLOTYPECALLER as GATK4_HAPLOTYPECALLER_PREPARE   } from '../../../modules/nf-core/gatk4/haplotypecaller/main'
include { GATK4_SELECTVARIANTS as GATK4_SELECTVARIANTS_SNP         } from '../../../modules/nf-core/gatk4/selectvariants/main'
include { GATK4_SELECTVARIANTS as GATK4_SELECTVARIANTS_INDEL       } from '../../../modules/nf-core/gatk4/selectvariants/main'
include { GATK4_VARIANTFILTRATION as GATK4_VARIANTFILTRATION_SNP   } from '../../../modules/nf-core/gatk4/variantfiltration/main'
include { GATK4_VARIANTFILTRATION as GATK4_VARIANTFILTRATION_INDEL } from '../../../modules/nf-core/gatk4/variantfiltration/main'
include { GATK4_SELECTVARIANTS as GATK4_SELECTVARIANTS_SNP_BQ      } from '../../../modules/nf-core/gatk4/selectvariants/main'
include { GATK4_SELECTVARIANTS as GATK4_SELECTVARIANTS_INDEL_BQ    } from '../../../modules/nf-core/gatk4/selectvariants/main'
include { GATK4_BASERECALIBRATOR                                   } from '../../../modules/nf-core/gatk4/baserecalibrator/main'
include { GATK4_APPLYBQSR                                          } from '../../../modules/nf-core/gatk4/applybqsr/main'
include { GATK4_BASERECALIBRATOR as GATK4_BASERECALIBRATOR_FORSTAT } from '../../../modules/nf-core/gatk4/baserecalibrator/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_BQSR                    } from '../../../modules/nf-core/samtools/index/main'
include { GATK4_ANALYZECOVARIATES                                  } from '../../../modules/local/gatk4_analyzecovariates'
include { GATK4_HAPLOTYPECALLER as GATK4_HAPLOTYPECALLER_GVCF      } from '../../../modules/nf-core/gatk4/haplotypecaller/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_FIL20                     } from '../../../modules/nf-core/samtools/view/main'
include { PICARD_SORTSAM                                           } from '../../../modules/nf-core/picard/sortsam/main'
include { SAMTOOLS_FAIDX                                           } from '../../../modules/nf-core/samtools/faidx/main'
include { PICARD_MARKDUPLICATES                                    } from '../../../modules/nf-core/picard/markduplicates/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_VARIANT_CALLING         } from '../../../modules/nf-core/samtools/index/main'
include { PICARD_COLLECTMULTIPLEMETRICS                            } from '../../../modules/nf-core/picard/collectmultiplemetrics/main'
include { PICARD_COLLECTINSERTSIZEMETRICS                          } from '../../../modules/nf-core/picard/collectinsertsizemetrics/main'
include { SAMTOOLS_DEPTH                                           } from '../../../modules/nf-core/samtools/depth/main'
include { PICARD_ADDORREPLACEREADGROUPS                            } from '../../../modules/nf-core/picard/addorreplacereadgroups/main' 
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_MAPPING                 } from '../../../modules/nf-core/samtools/index/main'
include { QUALIMAP_BAMQC                                           } from '../../../modules/nf-core/qualimap/bamqc/main'

include { GET_CHROMSOME_LST                                               } from '../../../modules/local/get_chromsomelst'
include { GATK4_GENOMICSDBIMPORT                                          } from '../../../modules/nf-core/gatk4/genomicsdbimport/main'
include { GATK4_GENOTYPEGVCFS                                             } from '../../../modules/nf-core/gatk4/genotypegvcfs/main'
include { GATK4_GATHERVCFS                                                } from '../../../modules/local/gatk4_gathervcfs'
include { GATK4_SELECTVARIANTS as GATK4_SELECTVARIANTS_SNP_HQ             } from '../../../modules/nf-core/gatk4/selectvariants/main'
include { GATK4_SELECTVARIANTS as GATK4_SELECTVARIANTS_INDEL_HQ           } from '../../../modules/nf-core/gatk4/selectvariants/main'
include { GATK4_VARIANTFILTRATION as GATK4_VARIANTFILTRATION_SNP_FOR_VQ   } from '../../../modules/nf-core/gatk4/variantfiltration/main'
include { GATK4_VARIANTFILTRATION as GATK4_VARIANTFILTRATION_INDEL_FOR_VQ } from '../../../modules/nf-core/gatk4/variantfiltration/main'
include { FILTER_PASS_VARIANT as FILTER_PASS_VARIANT_FOR_VQ               } from '../../../modules/local/filter_pass_variant'
include { GATK4_INDEXFEATUREFILE as GATK4_INDEXFEATUREFILE_FOR_VQ         } from '../../../modules/nf-core/gatk4/indexfeaturefile/main'
include { GATK4_VARIANTRECALIBRATOR as GATK4_VARIANTRECALIBRATOR_SNP      } from '../../../modules/nf-core/gatk4/variantrecalibrator/main'
include { GATK4_APPLYVQSR as GATK4_APPLYVQSR_SNP                          } from '../../../modules/nf-core/gatk4/applyvqsr/main'
include { GATK4_VARIANTRECALIBRATOR as GATK4_VARIANTRECALIBRATOR_INDEL    } from '../../../modules/nf-core/gatk4/variantrecalibrator/main'
include { GATK4_APPLYVQSR as GATK4_APPLYVQSR_INDEL                        } from '../../../modules/nf-core/gatk4/applyvqsr/main'
include { GATK4_SELECTVARIANTS as GATK4_SELECTVARIANTS_SNP_VQ             } from '../../../modules/nf-core/gatk4/selectvariants/main'
include { GATK4_SELECTVARIANTS as GATK4_SELECTVARIANTS_INDEL_VQ           } from '../../../modules/nf-core/gatk4/selectvariants/main'
include { FILTER_PASS_VARIANT as FILTER_PASS_VARIANT_VQ                   } from '../../../modules/local/filter_pass_variant'
include { GATK4_INDEXFEATUREFILE_GCVF                                     } from '../../../modules/nf-core/gatk4/indexgcvffile/main'









workflow VARIANT{

        take:
        ch_bwamem2_men_bam
	ch_fasta       

        main:
        def bwaidx_meta='fasta'
	ch_fasta_genome_withmeta=Channel.of(bwaidx_meta).combine(ch_fasta).collect().map{it}
        ch_fasta_genome_withoutmeta=ch_fasta.view()

        //ch_fasta_genome_withmeta=tuple(bwaidx_meta, file(params.fasta))
        //ch_fasta_genome_withoutmeta=file(params.fasta)

        PICARD_ADDORREPLACEREADGROUPS(
        ch_bwamem2_men_bam
        )
        ch_bam_to_fil = PICARD_ADDORREPLACEREADGROUPS.out.bam

        SAMTOOLS_INDEX_MAPPING(
            PICARD_ADDORREPLACEREADGROUPS.out.bam
        )
        ch_samtools_index_mapping_bai = SAMTOOLS_INDEX_MAPPING.out.bai


        QUALIMAP_BAMQC(
            PICARD_ADDORREPLACEREADGROUPS.out.bam,
            file("NOFILE")
        )


            // start sample variants calling
        ch_bam_to_fil=PICARD_ADDORREPLACEREADGROUPS.out.bam
        ch_bam_bai_to_fil=ch_bam_to_fil.join(ch_samtools_index_mapping_bai)       

        SAMTOOLS_VIEW_FIL20(
            ch_bam_bai_to_fil,
            ch_fasta_genome_withmeta,
            file("NOFILE")
        )

        PICARD_SORTSAM(
            SAMTOOLS_VIEW_FIL20.out.bam,
            "coordinate"
        )
        
        SAMTOOLS_FAIDX(
            ch_fasta_genome_withmeta
        )
        ch_fai_genome_withoutmeta=SAMTOOLS_FAIDX.out.fai.collect().flatten().toList().map{it -> it[1]}

        PICARD_MARKDUPLICATES(
            PICARD_SORTSAM.out.bam,
            ch_fasta_genome_withmeta,
            SAMTOOLS_FAIDX.out.fai
        )

        SAMTOOLS_INDEX_VARIANT_CALLING(
            PICARD_MARKDUPLICATES.out.bam
        )
        
        ch_bam_to_calling=PICARD_MARKDUPLICATES.out.bam
        ch_bam_bai_to_calling=ch_bam_to_calling.join(SAMTOOLS_INDEX_VARIANT_CALLING.out.bai)
        PICARD_COLLECTMULTIPLEMETRICS(
            ch_bam_bai_to_calling,
            ch_fasta_genome_withmeta,
            SAMTOOLS_FAIDX.out.fai        
        )

        PICARD_COLLECTINSERTSIZEMETRICS(
            ch_bam_to_calling       
        )

        depth_meta2="intervals"
        SAMTOOLS_DEPTH(
            ch_bam_to_calling,
            tuple(depth_meta2, file("NOFILE"))
        )
        GATK4_CREATESEQUENCEDICTIONARY(
            ch_fasta_genome_withmeta
        )
        ch_dict_genome_withoutmeta=GATK4_CREATESEQUENCEDICTIONARY.out.dict.collect().flatten().toList().map{it -> it[1]}

        //haplotypecaller_intervals=tuple(file("NOFILE"),file("NOFILE"))
        ch_hap_prepare_calling=ch_bam_bai_to_calling.map{it -> tuple(it[0],it[1],it[2],file("NOFILE1"),file("NOFILE2"))}
        GATK4_HAPLOTYPECALLER_PREPARE(
            ch_hap_prepare_calling,
            ch_fasta_genome_withoutmeta,
            ch_fai_genome_withoutmeta,
            ch_dict_genome_withoutmeta,
            file("NOFILE4"),
            file("NOFILE5")
        )
        ch_vcf_tbi_selectvariants=GATK4_HAPLOTYPECALLER_PREPARE.out.vcf.join(GATK4_HAPLOTYPECALLER_PREPARE.out.tbi)
        ch_vcf_tbi_intervals_selectvariants=ch_vcf_tbi_selectvariants.map{
            it -> tuple(it[0],it[1],it[2],file("NOFILE"))
        }
        GATK4_SELECTVARIANTS_SNP(
            ch_vcf_tbi_intervals_selectvariants
        )

        GATK4_SELECTVARIANTS_INDEL(
            ch_vcf_tbi_intervals_selectvariants
        )
        
        ch_tbi_variantfil_snp=GATK4_SELECTVARIANTS_SNP.out.vcf
        ch_vcf_tbi_variantfil_snp=ch_tbi_variantfil_snp.join(GATK4_SELECTVARIANTS_SNP.out.tbi)
        
        index_fai_variantfil='fai'
        index_dict_variantfil='dict'
        GATK4_VARIANTFILTRATION_SNP(
            ch_vcf_tbi_variantfil_snp,
            ch_fasta_genome_withmeta,
            SAMTOOLS_FAIDX.out.fai,
            GATK4_CREATESEQUENCEDICTIONARY.out.dict
        )

        ch_tbi_variantfil_indel=GATK4_SELECTVARIANTS_INDEL.out.vcf
        ch_vcf_tbi_variantfil_indel=ch_tbi_variantfil_indel.join(GATK4_SELECTVARIANTS_INDEL.out.tbi)
        GATK4_VARIANTFILTRATION_INDEL(
            ch_vcf_tbi_variantfil_indel,
            ch_fasta_genome_withmeta,
            SAMTOOLS_FAIDX.out.fai,
            GATK4_CREATESEQUENCEDICTIONARY.out.dict
        )

        ch_vcf_tbi_selectvariants_snp_forbq=GATK4_VARIANTFILTRATION_SNP.out.vcf.join(GATK4_VARIANTFILTRATION_SNP.out.tbi)
        ch_vcf_tbi_intervals_selectvariants_snp_forbq=ch_vcf_tbi_selectvariants_snp_forbq.map{
            it -> tuple(it[0],it[1],it[2],file("NOFILE"))
        }
        GATK4_SELECTVARIANTS_SNP_BQ(
            ch_vcf_tbi_intervals_selectvariants_snp_forbq
        )

        ch_vcf_tbi_selectvariants_indel_forbq=GATK4_VARIANTFILTRATION_INDEL.out.vcf.join(GATK4_VARIANTFILTRATION_INDEL.out.tbi)
        ch_vcf_tbi_intervals_selectvariants_indel_forbq=ch_vcf_tbi_selectvariants_indel_forbq.map{
            it -> tuple(it[0],it[1],it[2],file("NOFILE"))
        }
        GATK4_SELECTVARIANTS_INDEL_BQ(
            ch_vcf_tbi_intervals_selectvariants_indel_forbq
        )

        ch_vcf_bqsr=GATK4_SELECTVARIANTS_SNP_BQ.out.vcf.join(GATK4_SELECTVARIANTS_INDEL_BQ.out.vcf)
        ch_tbi_bqsr=GATK4_SELECTVARIANTS_SNP_BQ.out.tbi.join(GATK4_SELECTVARIANTS_INDEL_BQ.out.tbi)
        ch_bam_bqsr=PICARD_MARKDUPLICATES.out.bam.join(SAMTOOLS_INDEX_VARIANT_CALLING.out.bai)
        ch_bam_vcf_tbi_intervals_bqsr_join=ch_bam_bqsr.join(ch_vcf_bqsr).join(ch_tbi_bqsr)
        ch_bam_vcf_tbi_intervals_bqsr=ch_bam_vcf_tbi_intervals_bqsr_join.map{
                                            it -> 
                                            def vcfforbqsr = [it[3],it[4]]
                                            def tbiforbqsr = [it[5],it[6]]
                                            tuple(it[0],it[1],it[2],file("NOFILE1"),vcfforbqsr,tbiforbqsr)
                                        }
        GATK4_BASERECALIBRATOR(
            ch_bam_vcf_tbi_intervals_bqsr,
            ch_fasta_genome_withoutmeta,
            ch_fai_genome_withoutmeta,
            ch_dict_genome_withoutmeta
        )

        ch_bam_bai_table_intervals_applybqsr=ch_bam_bqsr.join(GATK4_BASERECALIBRATOR.out.table).map{it -> tuple(it[0],it[1],it[2],it[3],file("NOFILE1"))}
        GATK4_APPLYBQSR(
            ch_bam_bai_table_intervals_applybqsr,
            ch_fasta_genome_withoutmeta,
            ch_fai_genome_withoutmeta,
            ch_dict_genome_withoutmeta
        )
        SAMTOOLS_INDEX_BQSR(
            GATK4_APPLYBQSR.out.bam
        )

        ch_bam_bqsr_stat=GATK4_APPLYBQSR.out.bam.join(SAMTOOLS_INDEX_BQSR.out.bai)
        ch_bam_vcf_tbi_intervals_bqsrforstat_join=ch_bam_bqsr_stat.join(ch_vcf_bqsr).join(ch_tbi_bqsr)
        ch_bam_vcf_tbi_intervals_bqsrforstat=ch_bam_vcf_tbi_intervals_bqsrforstat_join.map{
                                            it -> 
                                            def vcfforbqsr = [it[3],it[4]]
                                            def tbiforbqsr = [it[5],it[6]]
                                            tuple(it[0],it[1],it[2],file("NOFILE1"),vcfforbqsr,tbiforbqsr)
                                        }
        GATK4_BASERECALIBRATOR_FORSTAT(
            ch_bam_vcf_tbi_intervals_bqsrforstat,
            ch_fasta_genome_withoutmeta,
            ch_fai_genome_withoutmeta,
            ch_dict_genome_withoutmeta
        )


        ch_before_after_bqsr_analysis=GATK4_BASERECALIBRATOR.out.table.join(GATK4_BASERECALIBRATOR_FORSTAT.out.table)
        GATK4_ANALYZECOVARIATES(
            ch_before_after_bqsr_analysis
        )

        haplotypecaller_intervals=tuple(file("NOFILE1"),file("NOFILE2"))
        ch_hap_gvcf_calling=ch_bam_bqsr_stat.map{it -> tuple(it[0],it[1],it[2],file("NOFILE1"),file("NOFILE2"))}
        GATK4_HAPLOTYPECALLER_GVCF(
            ch_hap_gvcf_calling,
            ch_fasta_genome_withoutmeta,
            ch_fai_genome_withoutmeta,
            ch_dict_genome_withoutmeta,
            file("NOFILE4"),
            file("NOFILE5")
        )
        //def bwaidx_meta='fasta'
        //ch_fasta_genome_withmeta=tuple(bwaidx_meta, file(params.fasta))
        //ch_fasta_genome_withoutmeta=file(params.fasta)

//        SAMTOOLS_FAIDX(
//                ch_fasta_genome_withmeta
//            )
//        ch_fai_genome_withoutmeta= SAMTOOLS_FAIDX.out.fai.collect().flatten().toList().map{it -> it[1]} 
        ch_samtools_faidx_fai    = SAMTOOLS_FAIDX.out.fai

//        GATK4_CREATESEQUENCEDICTIONARY(
//                ch_fasta_genome_withmeta
//            )
//        ch_dict_genome_withoutmeta=GATK4_CREATESEQUENCEDICTIONARY.out.dict.collect().flatten().toList().map{it -> it[1]}
        ch_gatk4_createsequencedictionary_dict = GATK4_CREATESEQUENCEDICTIONARY.out.dict

        GET_CHROMSOME_LST(
            ch_fasta_genome_withoutmeta
        )

        //ch_gvcf_indexfeaturefile = GATK4_HAPLOTYPECALLER_GVCF.out.vcf.map{ it -> it[1]}.collect().toList()
                           
        //GATK4_INDEXFEATUREFILE_GCVF(
        //    ch_gvcf_indexfeaturefile
        //)
        //ch_gvcf_dbimport.view() 
        //ch_tbi_dbimport = GATK4_INDEXFEATUREFILE_GCVF.out.tbi.collect().toList()
        // ch_tbi_dbimport.view()   
    
        
        //ch_gvcf_dbimport=GATK4_HAPLOTYPECALLER_GVCF.out.vcf.map{it -> it[1]}.collect{it -> "--variant ${it}"}.map{it-> it.join(' ')}
        ch_gvcf_dbimport=GATK4_HAPLOTYPECALLER_GVCF.out.vcf.map{it -> it[1]}.collect().toList()
        ch_tbi_dbimport=GATK4_HAPLOTYPECALLER_GVCF.out.tbi.map{it -> it[1]}.collect().toList()
        ch_gvcf_tbi_interval_wspace_dbimport=Channel.from(["id":"population"]).combine(
            ch_gvcf_dbimport).combine(ch_tbi_dbimport).combine(GET_CHROMSOME_LST.out.lst.map{it->it[0]}).combine(
            Channel.from(false)).combine(Channel.from(file("NOFILE")))
        GATK4_GENOMICSDBIMPORT(
            ch_gvcf_tbi_interval_wspace_dbimport,
            false,
            false,
            false
        )
    
        ch_chrom_genotypegvcfs=GET_CHROMSOME_LST.out.lst.map{it->it[1]}.splitText(by:1).map{it->it.split("\n")[0]}
        ch_ctg_genotypegvcfs=GET_CHROMSOME_LST.out.lst.map{it->it[2]}
        ch_chrom_ctg_genotypegvcfs=ch_chrom_genotypegvcfs.map{file("NOFILE")}
        ch_chrom_ctg_genotypegvcfs=ch_chrom_ctg_genotypegvcfs.concat(ch_ctg_genotypegvcfs)
        ch_ctgline_genotypegvcfs=Channel.from(false)
        ch_chrom_genotypegvcfs=ch_chrom_genotypegvcfs.concat(ch_ctgline_genotypegvcfs)
        ch_meta_id_genotypegvcfs=ch_chrom_genotypegvcfs.map{it->
            it==false ? tuple(["id":"ctg"],it) : tuple(["id":it],it)
        }
        ch_interval_genotypegvcfs=ch_meta_id_genotypegvcfs.merge(ch_chrom_ctg_genotypegvcfs)
        ch_gvcfdb_genotypegvcfs=GATK4_GENOMICSDBIMPORT.out.genomicsdb.map{itdb->itdb[1]}
        ch_gvcf_fai_interval_genotypegvcfs=ch_interval_genotypegvcfs.combine(ch_gvcf_dbimport).combine(ch_tbi_dbimport).combine(ch_gvcfdb_genotypegvcfs).map{
            it-> tuple(it[0],it[5],it[4],it[1],it[2])
        }
        //GATK4_GENOMICSDBIMPORT.out.genomicsdb.map{itdb->itdb[1]}.view()
    
        GATK4_GENOTYPEGVCFS(
            ch_gvcf_fai_interval_genotypegvcfs,
            ch_fasta_genome_withoutmeta,
            ch_fai_genome_withoutmeta,
            ch_dict_genome_withoutmeta
            //file("NOFILE"),
            //file("NOFILE")
        )
    
        ch_vcf_gathervcfs=ch_interval_genotypegvcfs.merge(
            Channel.from(1..100000)).join(
                GATK4_GENOTYPEGVCFS.out.vcf).toSortedList{
                     a, b -> a[3] <=> b[3] 
                }.flatten().buffer(size:5).map{it->it[4]}.flatten().collect().map{it-> tuple(['id':'population'],it)}
        ch_tbi_gathervcfs=ch_interval_genotypegvcfs.merge(
            Channel.from(1..100000)).join(
                GATK4_GENOTYPEGVCFS.out.tbi).toSortedList{
                     a, b -> a[3] <=> b[3] 
                }.flatten().buffer(size:5).map{it->it[4]}.flatten().collect().map{it-> tuple(['id':'population'],it)}
        ch_vcf_tbi_gathervcfs=ch_vcf_gathervcfs.join(ch_tbi_gathervcfs)
        GATK4_GATHERVCFS(
            ch_vcf_tbi_gathervcfs
        )

        ch_vcf_tbi_selectvariants_snp_hq=GATK4_GATHERVCFS.out.vcf.join(GATK4_GATHERVCFS.out.tbi)
        ch_vcf_tbi_intervals_selectvariants_snp_hq=ch_vcf_tbi_selectvariants_snp_hq.map{
                it -> tuple(it[0],it[1],it[2],file("NOFILE"))
        }

        GATK4_SELECTVARIANTS_SNP_HQ(
            ch_vcf_tbi_intervals_selectvariants_snp_hq
        )
        GATK4_SELECTVARIANTS_INDEL_HQ(
            ch_vcf_tbi_intervals_selectvariants_snp_hq
        )
       
        ch_vcf_tbi_variantfil_snp_forvq=GATK4_SELECTVARIANTS_SNP_HQ.out.vcf.join(GATK4_SELECTVARIANTS_SNP_HQ.out.tbi)
        GATK4_VARIANTFILTRATION_SNP_FOR_VQ(
            ch_vcf_tbi_variantfil_snp_forvq,
            ch_fasta_genome_withmeta,
            ch_samtools_faidx_fai,
            ch_gatk4_createsequencedictionary_dict
        )
        ch_vcf_tbi_variantfil_indel_forvq=GATK4_SELECTVARIANTS_INDEL_HQ.out.vcf.join(GATK4_SELECTVARIANTS_INDEL_HQ.out.tbi)
        GATK4_VARIANTFILTRATION_INDEL_FOR_VQ(
            ch_vcf_tbi_variantfil_indel_forvq,
            ch_fasta_genome_withmeta,
            ch_samtools_faidx_fai,
	    ch_gatk4_createsequencedictionary_dict   
        )

        ch_vcf_tbi_type_filpass_for_vq=GATK4_VARIANTFILTRATION_SNP_FOR_VQ.out.vcf.join(GATK4_VARIANTFILTRATION_SNP_FOR_VQ.out.tbi).map{
            it-> tuple(['id':'snp'],it[1],it[2])
        }.concat(GATK4_VARIANTFILTRATION_INDEL_FOR_VQ.out.vcf.join(GATK4_VARIANTFILTRATION_INDEL_FOR_VQ.out.tbi).map{
            it-> tuple(['id':'indel'],it[1],it[2])
        })
        FILTER_PASS_VARIANT_FOR_VQ(

            ch_vcf_tbi_type_filpass_for_vq
        )
        
        GATK4_INDEXFEATUREFILE_FOR_VQ(
            FILTER_PASS_VARIANT_FOR_VQ.out.vcf
        )
    
        ch_resource_vcf_vqsr_indel=FILTER_PASS_VARIANT_FOR_VQ.out.vcf.map{it->if(it[0].id=="indel"){it[1]}}.combine(GATK4_VARIANTFILTRATION_INDEL_FOR_VQ.out.vcf.map{it->it[1]})
        ch_resource_tbi_vqsr_indel=GATK4_INDEXFEATUREFILE_FOR_VQ.out.index.map{it->if(it[0].id=="indel"){it[1]}}.combine(GATK4_VARIANTFILTRATION_INDEL_FOR_VQ.out.tbi.map{it->it[1]})
        ch_vcf_tbi_vqsr_indel=GATK4_GATHERVCFS.out.vcf.join(GATK4_GATHERVCFS.out.tbi).map{it->tuple(['id':'indel'],it[1],it[2])}
        // ch_resource_tbi_vqsr_indel.view()
        GATK4_VARIANTRECALIBRATOR_INDEL(
            ch_vcf_tbi_vqsr_indel,
            ch_resource_vcf_vqsr_indel,
            ch_resource_tbi_vqsr_indel,
            ch_fasta_genome_withoutmeta,
            ch_fai_genome_withoutmeta,
            ch_dict_genome_withoutmeta
        )

        ch_applyvqsr_indel=ch_vcf_tbi_vqsr_indel.join(
            GATK4_VARIANTRECALIBRATOR_INDEL.out.recal).join(
            GATK4_VARIANTRECALIBRATOR_INDEL.out.idx).join(
            GATK4_VARIANTRECALIBRATOR_INDEL.out.tranches)
        GATK4_APPLYVQSR_INDEL(
            ch_applyvqsr_indel,
            ch_fasta_genome_withoutmeta,
            ch_fai_genome_withoutmeta,
            ch_dict_genome_withoutmeta
        )

        ch_resource_vcf_vqsr_snp=FILTER_PASS_VARIANT_FOR_VQ.out.vcf.map{it->if(it[0].id=="snp"){it[1]}}.combine(GATK4_VARIANTFILTRATION_SNP_FOR_VQ.out.vcf.map{it->it[1]})
        ch_resource_tbi_vqsr_snp=GATK4_INDEXFEATUREFILE_FOR_VQ.out.index.map{it->if(it[0].id=="snp"){it[1]}}.combine(GATK4_VARIANTFILTRATION_SNP_FOR_VQ.out.tbi.map{it->it[1]})    
        ch_vcf_tbi_vqsr_indel_snp=GATK4_APPLYVQSR_INDEL.out.vcf.join(GATK4_APPLYVQSR_INDEL.out.tbi).map{it->tuple(['id':'indel_snp'],it[1],it[2])}
        GATK4_VARIANTRECALIBRATOR_SNP(
            ch_vcf_tbi_vqsr_indel_snp,
            ch_resource_vcf_vqsr_snp,
            ch_resource_tbi_vqsr_snp,
            ch_fasta_genome_withoutmeta,
            ch_fai_genome_withoutmeta,
            ch_dict_genome_withoutmeta

        )

        ch_applyvqsr_indel_snp=ch_vcf_tbi_vqsr_indel_snp.join(
            GATK4_VARIANTRECALIBRATOR_SNP.out.recal).join(
            GATK4_VARIANTRECALIBRATOR_SNP.out.idx).join(
            GATK4_VARIANTRECALIBRATOR_SNP.out.tranches)    
        GATK4_APPLYVQSR_SNP(
            ch_applyvqsr_indel_snp,
            ch_fasta_genome_withoutmeta,
            ch_fai_genome_withoutmeta,
            ch_dict_genome_withoutmeta
        )

        ch_vcf_tbi_intervals_selectvariants_snp_vq=GATK4_APPLYVQSR_SNP.out.vcf.join(GATK4_APPLYVQSR_SNP.out.tbi).map{
            it-> tuple(['id':'snp'],it[1],it[2],file("NOFILE"))
        }
        ch_vcf_tbi_intervals_selectvariants_indel_vq=GATK4_APPLYVQSR_SNP.out.vcf.join(GATK4_APPLYVQSR_SNP.out.tbi).map{
            it-> tuple(['id':'indel'],it[1],it[2],file("NOFILE"))
        }
        GATK4_SELECTVARIANTS_SNP_VQ(
            ch_vcf_tbi_intervals_selectvariants_snp_vq
        )
        GATK4_SELECTVARIANTS_INDEL_VQ(
            ch_vcf_tbi_intervals_selectvariants_indel_vq
        )

        ch_vcf_tbi_fil_pass_variant_vq=GATK4_SELECTVARIANTS_SNP_VQ.out.vcf.join(GATK4_SELECTVARIANTS_SNP_VQ.out.tbi).concat(
            GATK4_SELECTVARIANTS_INDEL_VQ.out.vcf.join(GATK4_SELECTVARIANTS_INDEL_VQ.out.tbi)).concat(
            GATK4_APPLYVQSR_SNP.out.vcf.join(GATK4_APPLYVQSR_SNP.out.tbi))
        FILTER_PASS_VARIANT_VQ(
            ch_vcf_tbi_fil_pass_variant_vq
        )

    emit:
       ch_vcf_gatk_indexfeaturefile_vq = FILTER_PASS_VARIANT_VQ.out.vcf
}
