/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def valid_params = [
    aligners       : ['star_salmon', 'star_rsem', 'hisat2'],
    trimmers       : ['trimgalore', 'fastp'],
    pseudoaligners : ['salmon'],
    rseqc_modules  : ['bam_stat', 'inner_distance', 'infer_experiment', 'junction_annotation', 'junction_saturation', 'read_distribution', 'read_duplication', 'tin']
]

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters

WorkflowRnaseq.initialise(params, log, valid_params)

// Check input path parameters to see if they exist
/*
checkPathParamList = [
    params.input, params.multiqc_config,
    params.fasta, params.transcript_fasta, params.additional_fasta,
    params.gtf, params.gff, params.gene_bed,
    params.ribo_database_manifest, params.splicesites,
    params.star_index, params.hisat2_index, params.rsem_index, params.salmon_index
]
*/
//for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
// if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// Check rRNA databases for sortmerna
if (params.remove_ribo_rna) {
    ch_ribo_db = file(params.ribo_database_manifest, checkIfExists: true)
    if (ch_ribo_db.isEmpty()) {exit 1, "File provided with --ribo_database_manifest is empty: ${ch_ribo_db.getName()}!"}
}

// Check if file with list of fastas is provided when running BBSplit
if (!params.skip_bbsplit && !params.bbsplit_index && params.bbsplit_fasta_list) {
    ch_bbsplit_fasta_list = file(params.bbsplit_fasta_list, checkIfExists: true)
    if (ch_bbsplit_fasta_list.isEmpty()) {exit 1, "File provided with --bbsplit_fasta_list is empty: ${ch_bbsplit_fasta_list.getName()}!"}
}

// Check alignment parameters
def prepareToolIndices  = []
if (!params.skip_bbsplit) { prepareToolIndices << 'bbsplit' }
if (!params.skip_alignment) { prepareToolIndices << params.aligner }
if (!params.skip_pseudo_alignment && params.pseudo_aligner) { prepareToolIndices << params.pseudo_aligner }

// Get RSeqC modules to run
def rseqc_modules = params.rseqc_modules ? params.rseqc_modules.split(',').collect{ it.trim().toLowerCase() } : []
if (params.bam_csi_index) {
    for (rseqc_module in ['read_distribution', 'inner_distance', 'tin']) {
        if (rseqc_modules.contains(rseqc_module)) {
            rseqc_modules.remove(rseqc_module)
        }
    }
}

// Stage dummy file to be used as an optional input where required
ch_dummy_file = file("$projectDir/assets/dummy_file.txt", checkIfExists: true)

// Check if an AWS iGenome has been provided to use the appropriate version of STAR
def is_aws_igenome = false
if (params.fasta && params.gtf) {
    if ((file(params.fasta).getName() - '.gz' == 'genome.fa') && (file(params.gtf).getName() - '.gz' == 'genes.gtf')) {
        is_aws_igenome = true
    }
}

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

// Header files for MultiQC
ch_pca_header_multiqc        = file("$projectDir/assets/multiqc/deseq2_pca_header.txt", checkIfExists: true)
ch_clustering_header_multiqc = file("$projectDir/assets/multiqc/deseq2_clustering_header.txt", checkIfExists: true)
ch_biotypes_header_multiqc   = file("$projectDir/assets/multiqc/biotypes_header.txt", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//
include { DUPRADAR                           } from '../modules/local/dupradar'
include { MULTIQC                            } from '../modules/local/multiqc'
include { MULTIQC_CUSTOM_BIOTYPE             } from '../modules/local/multiqc_custom_biotype'
include { UMITOOLS_PREPAREFORRSEM as UMITOOLS_PREPAREFORSALMON } from '../modules/local/umitools_prepareforrsem.nf'
include { FEATURECOUNTS                      } from '../modules/local/featurecounts'
include { DESEQ2_COUNTS                      } from '../modules/local/deseq2_featurecounts'
include { MAKE_GENE_INFO                     } from '../modules/local/make_gene_info'


//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { FASTQC_WORKFLOW         } from '../subworkflows/nf-core-test/FASTQC/main'
include { ALIGN          } from '../subworkflows/nf-core-test/ALIGN/main'
include { QUANTIFY       } from '../subworkflows/nf-core-test/quantity/main'
include { DIFFERENCE_EXPRESSION } from '../subworkflows/nf-core-test/DEG/main'
include { TWAS           } from '../subworkflows/nf-core-test/TWAS/main'


include { INPUT_CHECK    } from '../subworkflows/local/input_check'
include { PREPARE_GENOME } from '../subworkflows/local/prepare_genome'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { BBMAP_BBSPLIT               } from '../modules/nf-core/bbmap/bbsplit/main'
include { SUBREAD_FEATURECOUNTS       } from '../modules/nf-core/subread/featurecounts/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

//
// SUBWORKFLOW: Consisting entirely of nf-core/modules
//


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report     = []
def pass_mapped_reads  = [:]
def pass_trimmed_reads = [:]
def pass_strand_check  = [:]

workflow RNASEQ {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Uncompress and prepare reference genome files
    //
    def biotype = params.gencode ? "gene_type" : params.gtf_group_features
    PREPARE_GENOME (
        params.fasta,
        params.gtf,
        params.gff,
        params.splicesites,
        params.bbsplit_fasta_list,
        params.star_index,
        params.hisat2_index,
        params.bbsplit_index,
        is_aws_igenome,
        prepareToolIndices
    )
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)
    

    //
    // SUBWORKFLOW: Sub-sample FastQ files and pseudo-align with Salmon to auto-infer strandedness
    //
    // Return empty channel if ch_strand_fastq.auto_strand is empty so salmon index isn't created

    ch_fasta = PREPARE_GENOME.out.fasta
	
    // Check if contigs in genome fasta file > 512 Mbp
    if (!params.skip_alignment && !params.bam_csi_index) {
        PREPARE_GENOME
            .out
            .fai
            .map { WorkflowRnaseq.checkMaxContigSize(it, log) }
    }

    // Obtain input analysis process

    processList = params.processList.split(',')

    start_process = processList[0] //Obtain start process 

    if(start_process == "A"){
	if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
        FASTQC_WORKFLOW(
            ch_input,
            ch_fasta,
            PREPARE_GENOME.out.gtf,
            PREPARE_GENOME.out.bbsplit_index,
     	    prepareToolIndices,
            pass_trimmed_reads
        )
	ch_filtered_reads = FASTQC_WORKFLOW.out.ch_filtered_reads
        ch_versions = FASTQC_WORKFLOW.out.ch_versions
	
    }
    if(start_process == "B"){
        Channel.fromPath(params.input_algin)
        .splitCsv(header:true,sep:',')
        .map{
            create_fastq_channel(it)
        }
        .set{ch_filtered_reads}
	ALIGN(
            ch_filtered_reads,
            PREPARE_GENOME.out.star_index,
            PREPARE_GENOME.out.gtf,
            PREPARE_GENOME.out.fasta.map { [ [:], it ] },
	    PREPARE_GENOME.out.hisat2_index.map { [ [:], it ] },
            PREPARE_GENOME.out.splicesites.map { [ [:], it ] },
            ch_versions,
            is_aws_igenome
        )
    }else if("B" in processList ){
        
        ALIGN(
            ch_filtered_reads,
            PREPARE_GENOME.out.star_index,
            PREPARE_GENOME.out.gtf,
            PREPARE_GENOME.out.fasta.map { [ [:], it ] },
	    PREPARE_GENOME.out.hisat2_index.map { [ [:], it ] },
            PREPARE_GENOME.out.splicesites.map { [ [:], it ] },
	    ch_versions,
            is_aws_igenome
        )    
    }
    if(start_process == "C"){// 开始定量，定量分为转录本定量和基因定量
	Channel.fromPath(params.input_bam)
        .splitCsv(header:true,sep:',')
        .map{
            row -> create_bam_channel(row)
        }
	.set{ch_genome_bam}
	
	ch_genome_bam
	.map{row -> row[1]}.collect()
	.set{ch_genome_bam_total}
       
	
        QUANTIFY(
            PREPARE_GENOME.out.gtf,
            ch_genome_bam_total,
            ch_genome_bam,
	    biotype
        )
        ch_deseq_featurecounts = QUANTIFY.out.ch_deseq_featurecounts
        ch_tpm_pca            = QUANTIFY.out.ch_tpm_pca
        ch_tpm                = QUANTIFY.out.ch_tpm
    }else if("C" in processList ){ // 开始定量，定量分为转录本定量和基因定量
         
        QUANTIFY(
            PREPARE_GENOME.out.gtf,
            ALIGN.out.ch_genome_bam_total,
            ALIGN.out.ch_genome_bam,
	    biotype
        ) 
        ch_deseq_featurecounts = QUANTIFY.out.ch_deseq_featurecounts
        ch_tpm_pca             = QUANTIFY.out.ch_tpm_pca 
        ch_tpm                 = QUANTIFY.out.ch_tpm  
    }
    
    if(start_process == "D"){
	ch_tpm_pca		     = [[group:"counts"],[params.input_counts],[params.input_group]]
        ch_deseq_featurecounts   = [[group:"counts"],[params.input_counts],[params.input_compare],[params.input_group]]
       	 DIFFERENCE_EXPRESSION(
            	ch_deseq_featurecounts,
            	ch_tpm_pca
        )
    }else if("D" in processList){ // 开始进行差异表达分析
        
        ch_deseq_featurecounts.view()
        ch_tpm_pca.view() 
        DIFFERENCE_EXPRESSION(
            ch_deseq_featurecounts,
            ch_tpm_pca
            
        )
    }
    if(start_process == "E"){
        test_param  =  channel.of([["id":"twas"],[params.input_expression],[params.input_GM],[params.input_phenotype]])
        TWAS(
         test_param 
        )
    }else if("E" in processList){
       twas_param  =  ch_tpm.combine(Channel.of(params.pheno)).map{it -> tuple(["id":"twas"],it[1],it[2],it[3])}
       
        TWAS(
           twas_param
        )
    }	
     if("F" in processList){
     //
    // MODULE: Genome-wide coverage with BEDTools
    //
    if (!params.skip_alignment && !params.skip_bigwig) {
	
        BEDTOOLS_GENOMECOV (
            ch_genome_bam
        )
        ch_versions = ch_versions.mix(BEDTOOLS_GENOMECOV.out.versions.first())

        //
        // SUBWORKFLOW: Convert bedGraph to bigWig
        //
        BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_FORWARD (
            BEDTOOLS_GENOMECOV.out.bedgraph_forward,
            PREPARE_GENOME.out.chrom_sizes
        )
        ch_versions = ch_versions.mix(BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_FORWARD.out.versions)

        BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_REVERSE (
            BEDTOOLS_GENOMECOV.out.bedgraph_reverse,
            PREPARE_GENOME.out.chrom_sizes
        )
    }
   
    //
    // MODULE: Downstream QC steps
    //
    ch_dupradar_multiqc           = Channel.empty()
    ch_bamstat_multiqc            = Channel.empty()
    ch_inferexperiment_multiqc    = Channel.empty()
    ch_innerdistance_multiqc      = Channel.empty()
    ch_junctionannotation_multiqc = Channel.empty()
    ch_junctionsaturation_multiqc = Channel.empty()
    ch_readdistribution_multiqc   = Channel.empty()
    ch_readduplication_multiqc    = Channel.empty()
    ch_fail_strand_multiqc        = Channel.empty()
    ch_tin_multiqc                = Channel.empty()
    if (!params.skip_alignment && !params.skip_qc) {

        if (!params.skip_dupradar) {
            DUPRADAR (
                ch_genome_bam,
                PREPARE_GENOME.out.gtf
            )
            ch_dupradar_multiqc = DUPRADAR.out.multiqc
            ch_versions = ch_versions.mix(DUPRADAR.out.versions.first())
        }

        if (!params.skip_rseqc && rseqc_modules.size() > 0) {
            BAM_RSEQC (
                ch_genome_bam.join(ch_genome_bam_index, by: [0]),
                PREPARE_GENOME.out.gene_bed,
                rseqc_modules
            )
            ch_bamstat_multiqc            = BAM_RSEQC.out.bamstat_txt
            ch_inferexperiment_multiqc    = BAM_RSEQC.out.inferexperiment_txt
            ch_innerdistance_multiqc      = BAM_RSEQC.out.innerdistance_freq
            ch_junctionannotation_multiqc = BAM_RSEQC.out.junctionannotation_log
            ch_junctionsaturation_multiqc = BAM_RSEQC.out.junctionsaturation_rscript
            ch_readdistribution_multiqc   = BAM_RSEQC.out.readdistribution_txt
            ch_readduplication_multiqc    = BAM_RSEQC.out.readduplication_pos_xls
            ch_tin_multiqc                = BAM_RSEQC.out.tin_txt
            ch_versions = ch_versions.mix(BAM_RSEQC.out.versions)

            ch_inferexperiment_multiqc
                .map {
                    meta, strand_log ->
                        def inferred_strand = WorkflowRnaseq.getInferexperimentStrandedness(strand_log, 30)
                        pass_strand_check[meta.id] = true
                        if (meta.strandedness != inferred_strand[0]) {
                            pass_strand_check[meta.id] = false
                            return [ "$meta.id\t$meta.strandedness\t${inferred_strand.join('\t')}" ]
                        }
                }
                .collect()
                .map {
                    tsv_data ->
                        def header = [
                            "Sample",
                            "Provided strandedness",
                            "Inferred strandedness",
                            "Sense (%)",
                            "Antisense (%)",
                            "Undetermined (%)"
                        ]
                        WorkflowRnaseq.multiqcTsvFromList(tsv_data, header)
                }
                .set { ch_fail_strand_multiqc }
        }
    }

    //
    // SUBWORKFLOW: Pseudo-alignment and quantification with Salmon
    //
    ch_salmon_multiqc                   = Channel.empty()
    ch_pseudoaligner_pca_multiqc        = Channel.empty()
    ch_pseudoaligner_clustering_multiqc = Channel.empty()
    if (!params.skip_pseudo_alignment && params.pseudo_aligner == 'salmon') {
        QUANTIFY_SALMON (
            ch_filtered_reads,
            PREPARE_GENOME.out.salmon_index,
            ch_dummy_file,
            PREPARE_GENOME.out.gtf,
            false,
            params.salmon_quant_libtype ?: ''
        )
        ch_salmon_multiqc = QUANTIFY_SALMON.out.results
        ch_versions = ch_versions.mix(QUANTIFY_SALMON.out.versions)

        if (!params.skip_qc & !params.skip_deseq2_qc) {
            DESEQ2_QC_SALMON (
                QUANTIFY_SALMON.out.counts_gene_length_scaled,
                ch_pca_header_multiqc,
                ch_clustering_header_multiqc
            )
            ch_pseudoaligner_pca_multiqc        = DESEQ2_QC_SALMON.out.pca_multiqc
            ch_pseudoaligner_clustering_multiqc = DESEQ2_QC_SALMON.out.dists_multiqc
            ch_versions = ch_versions.mix(DESEQ2_QC_SALMON.out.versions)
        }
    }

    //
    // MODULE: Pipeline reporting
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    if (!params.skip_multiqc) {
        workflow_summary    = WorkflowRnaseq.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        methods_description    = WorkflowRnaseq.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
        ch_methods_description = Channel.value(methods_description)

        MULTIQC (
            ch_multiqc_config,
            ch_multiqc_custom_config.collect().ifEmpty([]),
            CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect(),
            ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
            ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'),
            ch_multiqc_logo.collect().ifEmpty([]),
            ch_fail_trimming_multiqc.collectFile(name: 'fail_trimmed_samples_mqc.tsv').ifEmpty([]),
            ch_fail_mapping_multiqc.collectFile(name: 'fail_mapped_samples_mqc.tsv').ifEmpty([]),
            ch_fail_strand_multiqc.collectFile(name: 'fail_strand_check_mqc.tsv').ifEmpty([]),
            ch_fastqc_raw_multiqc.collect{it[1]}.ifEmpty([]),
            ch_fastqc_trim_multiqc.collect{it[1]}.ifEmpty([]),
            ch_trim_log_multiqc.collect{it[1]}.ifEmpty([]),
            ch_sortmerna_multiqc.collect{it[1]}.ifEmpty([]),
            ch_star_multiqc.collect{it[1]}.ifEmpty([]),
            ch_hisat2_multiqc.collect{it[1]}.ifEmpty([]),
            ch_rsem_multiqc.collect{it[1]}.ifEmpty([]),
            ch_salmon_multiqc.collect{it[1]}.ifEmpty([]),
            ch_samtools_stats.collect{it[1]}.ifEmpty([]),
            ch_samtools_flagstat.collect{it[1]}.ifEmpty([]),
            ch_samtools_idxstats.collect{it[1]}.ifEmpty([]),
            ch_markduplicates_multiqc.collect{it[1]}.ifEmpty([]),
            ch_featurecounts_multiqc.collect{it[1]}.ifEmpty([]),
            ch_aligner_pca_multiqc.collect().ifEmpty([]),
            ch_aligner_clustering_multiqc.collect().ifEmpty([]),
            ch_pseudoaligner_pca_multiqc.collect().ifEmpty([]),
            ch_pseudoaligner_clustering_multiqc.collect().ifEmpty([]),
            ch_preseq_multiqc.collect{it[1]}.ifEmpty([]),
            ch_qualimap_multiqc.collect{it[1]}.ifEmpty([]),
            ch_dupradar_multiqc.collect{it[1]}.ifEmpty([]),
            ch_bamstat_multiqc.collect{it[1]}.ifEmpty([]),
            ch_inferexperiment_multiqc.collect{it[1]}.ifEmpty([]),
            ch_innerdistance_multiqc.collect{it[1]}.ifEmpty([]),
            ch_junctionannotation_multiqc.collect{it[1]}.ifEmpty([]),
            ch_junctionsaturation_multiqc.collect{it[1]}.ifEmpty([]),
            ch_readdistribution_multiqc.collect{it[1]}.ifEmpty([]),
            ch_readduplication_multiqc.collect{it[1]}.ifEmpty([]),
            ch_tin_multiqc.collect{it[1]}.ifEmpty([])
        )
        multiqc_report = MULTIQC.out.report.toList()
    }
  }
	// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
    def create_fastq_channel(LinkedHashMap row) {
        // create meta map
        def meta = [:]
        meta.id           = row.sample
        

        // add path(s) of the fastq file(s) to the meta map
        def fastq_meta = []
        if (!file(row.fastq_1).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
        }else{
            meta.single_end  = true 
            meta.strandedness = row.strandedness
            fastq_meta = [ meta, [ file(row.fastq_1) ] ]
        }
        if (!file(row.fastq_2).exists()) {
                exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
            }else{
                meta.single_end  = false
                meta.strandedness = row.strandedness
                fastq_meta = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
            }

    return fastq_meta
    }

// Function to get list of [ meta, bam ]
    def create_bam_channel(LinkedHashMap row) {
        // create meta map
        def meta = [:]
        meta.id           = row.sample
        // add path(s) of the fastq file(s) to the meta map
        def fastq_meta = []
        if (row.single == true ) {
            meta.single_end  = true 
            meta.strandedness = row.strandedness
            fastq_meta = [ meta, file(row.bam)  ]
        }else{
                meta.single_end  = false 
                meta.strandedness = row.strandedness
                fastq_meta = [ meta, file(row.bam) ]
            }

    return fastq_meta
    }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report, pass_mapped_reads, pass_trimmed_reads, pass_strand_check)
    }

    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }

    NfcoreTemplate.summary(workflow, params, log, pass_mapped_reads, pass_trimmed_reads, pass_strand_check)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
