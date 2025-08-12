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
    ch_vcf
    //ch_tbi_dbimport
    main:
    def bwaidx_meta='fasta'
    ch_fasta_genome_withmeta=tuple(bwaidx_meta, file(params.fasta))
    ch_fasta_genome_withoutmeta=file(params.fasta)
    
    GATK4_INDEXFEATUREFILE_VQ(
        ch_vcf
    )


    emit:
    ch_vcf_tbi_annotable = ch_vcf.map{it-> if(it[0].id=='snp'){it}}.join(
        GATK4_INDEXFEATUREFILE_VQ.out.index.map{it-> if(it[0].id=='snp'){it}}).map{it-> tuple(['id':'ANNO'],it[1],it[2])}
    ch_vcf_tbi_structure_vcfaddid=ch_vcf.map{it-> if(it[0].id=='snp'){it}}.join( GATK4_INDEXFEATUREFILE_VQ.out.index.map{it-> if(it[0].id=='snp'){it}}).map{it-> tuple(['id':'popu_snp'],it[1],it[2])}  
    ch_vcf_tbi_pca_vcfaddid=ch_vcf.map{it-> if(it[0].id=='snp'){it}}.join(    GATK4_INDEXFEATUREFILE_VQ.out.index.map{it-> if(it[0].id=='snp'){it}}).map{it-> tuple(it[1],it[2])}
    ch_groups_vcf_emmax_groupping=ch_vcf.map{it-> if(it[0].id=='snp'){it}}.join(GATK4_INDEXFEATUREFILE_VQ.out.index.map{it-> if(it[0].id=='snp'){it}}).map{it-> tuple(['id':'popu_snp'],0,it[1],it[2])}
}
