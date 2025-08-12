/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FUSION_TWAS            } from '../subworkflows/local/fusion-twas/main'
include { EQTL                   } from '../subworkflows/local/eQTL/main'
include { EQTM                   } from '../subworkflows/local/eQTM/main'
include { MEQTL                  } from '../subworkflows/local/meQTL/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { PUREMETHLATION         } from '../subworkflows/local/pureMethlation/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow OMIC {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = Channel.empty()

    processList = params.processList.split(',')
    start_process = processList[0]

    
    if (start_process == "A"){
        EQTL()
    }
    if(start_process == "B"){
        FUSION_TWAS()
    }
    if(start_process == "C"){
        EQTM()
    }
    if(start_process == "D"){
        MEQTL()
    }
    if(start_process == "E"){
        PUREMETHLATION()
    }
    
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
