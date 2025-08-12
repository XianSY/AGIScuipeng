include { FASTQS     } from '../subworkflows/nf-core/fastqc'
include { ALIGN      } from '../subworkflows/nf-core/align'
include { VARIANT    } from '../subworkflows/nf-core/variant'
include { POPULATION } from '../subworkflows/nf-core/population'
include { ANNOTATION } from '../subworkflows/nf-core/annotation'
include { STRUCTURE  } from '../subworkflows/nf-core/structure'
include { PCA        } from '../subworkflows/nf-core/pca'
include { GWAS       } from '../subworkflows/nf-core/gwas'



include { GATK4_INDEXFEATUREFILE as GATK4_INDEXFEATUREFILE_VQ             } from '../modules/nf-core/gatk4/indexfeaturefile/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowTest.initialise(params, log)

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
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_MAPPING         } from '../modules/nf-core/samtools/index/main'




//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK                      } from '../subworkflows/local/input_check'
include { INPUT_CHECK as INPUT_CHECK_BAM   } from '../subworkflows/local/input_check'
include { PREPARE_GENOME                   } from '../subworkflows/local/prepare_genome' 
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
    def multiqc_report = []

workflow TEST {
	
    main:
    PREPARE_GENOME(
        params.fasta,
        params.gff,
        params.gtf,
        params.adapter_fasta
    )

    ch_versions = Channel.empty()
     
    processList = params.processList.split(',')
    start =  processList[0]

    if (start == "A"){
        INPUT_CHECK (
         params.input
   	     )
 	    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

        FASTQS(
                INPUT_CHECK.out.reads,
                params.save_trimmed_fail,
                params.save_merged
            )
        ch_merge_fastq_files = FASTQS.out.fastqc_clean.map{
        row -> def sampleParts = row[0].id.tokenize('_')
            tuple(sampleParts[0],sampleParts[1],row[0].single_end,row[1][0], row[1][1])
        }.groupTuple()

        ch_versions     = ch_versions.mix((FASTQS.out.versions))
        }
	if (start == "B"){
        //ch_merge_fastq_files = channel.fromPath(params.input).splitCsv().map{
        //row -> def sampleParts = row[0].tokenize('_')
        //    tuple(sampleParts[0],sampleParts[1],row[1],row[2],row[3])
        //}.groupTuple()
        ch_merge_fastq_files = channel.fromPath(params.input_align).splitCsv( header:true).map{
        row -> def sampleParts = row.sample.tokenize('_')
            tuple(sampleParts[0],sampleParts[1],row.single_end,row.fastq_1,row.fastq_2)
        }.groupTuple()
        ALIGN(
        ch_merge_fastq_files,
        PREPARE_GENOME.out.ch_fasta
        )
        ch_bwamem2_mem_bam          = ALIGN.out.bwamem2_mem_bam

        }else if('B' in processList){
            ALIGN(
                ch_merge_fastq_files,
                PREPARE_GENOME.out.ch_fasta
            )           
           ch_bwamem2_mem_bam          = ALIGN.out.bwamem2_mem_bam
        }
        if(start == "C"){
        ch_bwamem2_men_bam = channel.fromPath(params.input_bam)
                             .splitCsv(header:true)
                             .map{
                                row -> def bam=[:]
                                bam.id = row.sample
                                tuple(bam,row.bam)
                             }
        ch_bwamem2_men_bam.view()
        VARIANT(
            ch_bwamem2_men_bam          ,
            PREPARE_GENOME.out.ch_fasta
        )
        }else if('C' in processList){
            VARIANT(
                ch_bwamem2_mem_bam,
                PREPARE_GENOME.out.ch_fasta          
            )
        }
        if (start == 'D'){

	ch_vcf = Channel.fromPath(params.input_vcf).map{it -> tuple(['id':'snp'],it)}
	    GATK4_INDEXFEATUREFILE_VQ(ch_vcf) 

        ch_vcf_tbi_annotable = ch_vcf.map{it-> if(it[0].id=='snp'){it}}.join(
            GATK4_INDEXFEATUREFILE_VQ.out.index.map{it-> if(it[0].id=='snp'){it}}).map{it-> tuple(['id':'ANNO'],it[1],it[2])}
        ch_vcf_tbi_structure_vcfaddid=ch_vcf.map{it-> if(it[0].id=='snp'){it}}.join( GATK4_INDEXFEATUREFILE_VQ.out.index.map{it-> if(it[0].id=='snp'){it}}).map{it-> tuple(['id':'snp'],it[1],it[2])}  
        ch_vcf_tbi_pca_vcfaddid=ch_vcf.map{it-> if(it[0].id=='snp'){it}}.join(    GATK4_INDEXFEATUREFILE_VQ.out.index.map{it-> if(it[0].id=='snp'){it}}).map{it-> tuple(it[0],it[1],it[2])}
        ch_groups_vcf_emmax_groupping=ch_vcf.map{it-> if(it[0].id=='snp'){it}}.join(GATK4_INDEXFEATUREFILE_VQ.out.index.map{it-> if(it[0].id=='snp'){it}}).map{it-> tuple(['id':'snp'],0,it[1],it[2])}


        ANNOTATION(
            ch_vcf_tbi_annotable,
            PREPARE_GENOME.out.ch_gff,
            PREPARE_GENOME.out.ch_fasta
        )
        STRUCTURE(
           ch_vcf_tbi_structure_vcfaddid
        )
        PCA(
            ch_vcf_tbi_pca_vcfaddid,
            STRUCTURE.out.ch_plink_indep_prunein
            )
        GWAS(
            ch_vcf_tbi_pca_vcfaddid,
            ch_groups_vcf_emmax_groupping,
            ANNOTATION.out.ch_annoevf_emmax_anno
	        )
        }else if ('D' in processList){
            GATK4_INDEXFEATUREFILE_VQ(
            VARIANT.out.ch_vcf_gatk_indexfeaturefile_vq
        )

	ch_vcf_tbi_annotable = VARIANT.out.ch_vcf_gatk_indexfeaturefile_vq.map{it-> if(it[0].id=='snp'){it}}.join(
            GATK4_INDEXFEATUREFILE_VQ.out.index.map{it-> if(it[0].id=='snp'){it}}).map{it-> tuple(['id':'ANNO'],it[1],it[2])}
        ch_vcf_tbi_structure_vcfaddid=VARIANT.out.ch_vcf_gatk_indexfeaturefile_vq.map{it-> if(it[0].id=='snp'){it}}.join( GATK4_INDEXFEATUREFILE_VQ.out.index.map{it-> if(it[0].id=='snp'){it}}).map{it-> tuple(['id':'popu_snp'],it[1],it[2])}  
        ch_vcf_tbi_pca_vcfaddid=VARIANT.out.ch_vcf_gatk_indexfeaturefile_vq.map{it-> if(it[0].id=='snp'){it}}.join(    GATK4_INDEXFEATUREFILE_VQ.out.index.map{it-> if(it[0].id=='snp'){it}}).map{it-> tuple(['id':'popu_snp'],it[1],it[2])}
        ch_groups_vcf_emmax_groupping=VARIANT.out.ch_vcf_gatk_indexfeaturefile_vq.map{it-> if(it[0].id=='snp'){it}}.join(GATK4_INDEXFEATUREFILE_VQ.out.index.map{it-> if(it[0].id=='snp'){it}}).map{it-> tuple(['id':'popu_snp'],0,it[1],it[2])}



        ANNOTATION(
            ch_vcf_tbi_annotable,
            PREPARE_GENOME.out.ch_gff,
            PREPARE_GENOME.out.ch_fasta
        )
        STRUCTURE(
           ch_vcf_tbi_structure_vcfaddid
        )
       	PCA(
            ch_vcf_tbi_pca_vcfaddid,
        STRUCTURE.out.ch_plink_indep_prunein
        )
        GWAS(
          ch_vcf_tbi_pca_vcfaddid,
          ch_groups_vcf_emmax_groupping,
          ANNOTATION.out.ch_annoevf_emmax_anno
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
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log)
    }
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
