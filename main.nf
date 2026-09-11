#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/genomeqc
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/genomeqc
    Website: https://nf-co.re/genomeqc
    Slack  : https://nfcore.slack.com/channels/genomeqc
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GENOMEQC  } from './workflows/genomeqc'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_genomeqc_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_genomeqc_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_GENOMEQC {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    //
    // WORKFLOW: Run pipeline
    //
    GENOMEQC (
        samplesheet,
        params.multiqc_config,
        params.multiqc_logo,
        params.multiqc_methods_description,
        params.outdir,
        params.groups,
        params.busco_lineages_path,
        params.busco_lineage,
        params.gxdb,
        params.gxdb_manifest,
        params.ramdisk,
        params.repeat,
        params.skip_tidk,
        params.kvalue,
        params.RM_download_db,
        params.RM_db,
        params.famdb_library,
        params.famdb_lineage,
        params.run_repeatmodeler,
        params.te_clusterer,
        params.te,
        params.skip_busco,
        params.busco_config,
        params.busco_clean,
        params.ortho_version,
        params.val_tool,
        params.container_engine,
    )
    emit:
    multiqc_report = GENOMEQC.out.multiqc_report // channel: /path/to/multiqc_report.html
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_GENOMEQC (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        NFCORE_GENOMEQC.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
