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
include { GATK4_INDEXFEATUREFILE as GATK4_INDEXFEATUREFILE_VQ             } from '../../../modules/nf-core/gatk4/indexfeaturefile/main'
include { GATK4_CREATESEQUENCEDICTIONARY                           } from '../../../modules/nf-core/gatk4/createsequencedictionary/main'
include { SAMTOOLS_FAIDX                                           } from '../../../modules/nf-core/samtools/faidx/main'
include { GATK4_INDEXFEATUREFILE_GCVF                                     } from '../../../modules/nf-core/gatk4/indexgcvffile/main'






workflow POPULATION{
    take:
    ch_gvcf_dbimport
    //ch_tbi_dbimport
    main:
    def bwaidx_meta='fasta'
    ch_fasta_genome_withmeta=tuple(bwaidx_meta, file(params.fasta))
    ch_fasta_genome_withoutmeta=file(params.fasta)

    SAMTOOLS_FAIDX(
            ch_fasta_genome_withmeta
        )
    ch_fai_genome_withoutmeta= SAMTOOLS_FAIDX.out.fai.collect().flatten().toList().map{it -> it[1]} 
    ch_samtools_faidx_fai    = SAMTOOLS_FAIDX.out.fai

    GATK4_CREATESEQUENCEDICTIONARY(
            ch_fasta_genome_withmeta
        )
    ch_dict_genome_withoutmeta=GATK4_CREATESEQUENCEDICTIONARY.out.dict.collect().flatten().toList().map{it -> it[1]}
    ch_gatk4_createsequencedictionary_dict = GATK4_CREATESEQUENCEDICTIONARY.out.dict

    GET_CHROMSOME_LST(
        ch_fasta_genome_withoutmeta
    )
    ch_gvcf_indexfeaturefile = ch_gvcf_dbimport
                           .flatten()
    GATK4_INDEXFEATUREFILE_GCVF(
        ch_gvcf_indexfeaturefile
    )
    //ch_gvcf_dbimport.view() 
    ch_tbi_dbimport = GATK4_INDEXFEATUREFILE_GCVF.out.tbi.collect().toList()
    // ch_tbi_dbimport.view()   
    
    //ch_gvcf_dbimport=GATK4_HAPLOTYPECALLER_GVCF.out.vcf.map{it -> it[1]}.collect{it -> "--variant ${it}"}.map{it-> it.join(' ')}
    //ch_gvcf_dbimport=GATK4_HAPLOTYPECALLER_GVCF.out.vcf.map{it -> it[1]}.collect().toList()
    //ch_tbi_dbimport=GATK4_HAPLOTYPECALLER_GVCF.out.tbi.map{it -> it[1]}.collect().toList()
    ch_gvcf_tbi_interval_wspace_dbimport=Channel.from(["id":"population"]).combine(
        ch_gvcf_dbimport).combine(ch_tbi_dbimport).combine(GET_CHROMSOME_LST.out.lst.map{it->it[0]}).combine(
        Channel.from(false)).combine(Channel.from(file("NOFILE")))
    // ch_gvcf_tbi_interval_wspace_dbimport.view()
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
    //FILTER_PASS_VARIANT_FOR_VQ.out.vcf.view()
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
    
    GATK4_INDEXFEATUREFILE_VQ(
        FILTER_PASS_VARIANT_VQ.out.vcf
    )


    emit:
    ch_vcf_tbi_annotable = FILTER_PASS_VARIANT_VQ.out.vcf.map{it-> if(it[0].id=='snp'){it}}.join(
        GATK4_INDEXFEATUREFILE_VQ.out.index.map{it-> if(it[0].id=='snp'){it}}).map{it-> tuple(['id':'ANNO'],it[1],it[2])}
    ch_vcf_tbi_structure_vcfaddid=FILTER_PASS_VARIANT_VQ.out.vcf.map{it-> if(it[0].id=='snp'){it}}.join( GATK4_INDEXFEATUREFILE_VQ.out.index.map{it-> if(it[0].id=='snp'){it}}).map{it-> tuple(['id':'popu_snp'],it[1],it[2])}  
    ch_vcf_tbi_pca_vcfaddid=FILTER_PASS_VARIANT_VQ.out.vcf.map{it-> if(it[0].id=='snp'){it}}.join(    GATK4_INDEXFEATUREFILE_VQ.out.index.map{it-> if(it[0].id=='snp'){it}}).map{it-> tuple(it[1],it[2])}
    ch_groups_vcf_emmax_groupping=FILTER_PASS_VARIANT_VQ.out.vcf.map{it-> if(it[0].id=='snp'){it}}.join(GATK4_INDEXFEATUREFILE_VQ.out.index.map{it-> if(it[0].id=='snp'){it}}).map{it-> tuple(['id':'popu_snp'],0,it[1],it[2])}

}
