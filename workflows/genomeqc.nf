/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CREATE_PATH                       } from '../modules/local/create_path'
include { NCBIGENOMEDOWNLOAD                } from '../modules/nf-core/ncbigenomedownload/main'
include { BUSCO_BUSCO                       } from '../modules/nf-core/busco/busco/main'
include { GFFREAD                           } from '../modules/local/gffread'
include { ORTHOFINDER                       } from '../modules/nf-core/orthofinder/main'
include { FASTQC                            } from '../modules/nf-core/fastqc/main'
include { MULTIQC                           } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap                  } from 'plugin/nf-validation'
include { paramsSummaryMultiqc              } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML            } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText            } from '../subworkflows/local/utils_nfcore_genomeqc_pipeline'
include { validateInputSamplesheet          } from '../subworkflows/local/utils_nfcore_genomeqc_pipeline'
include { FASTA_EXPLORE_SEARCH_PLOT_TIDK    } from '../subworkflows/nf-core/fasta_explore_search_plot_tidk'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GENOMEQC {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    ch_samplesheet
        .map {
            validateInputSamplesheet(it) // Input validation
        }
        .branch {
            ncbi: it.size() == 2
            local: it.size() == 3
        }
        .set { ch_input }

    ch_input.ncbi.view()
    //
    // MODULE: Run create_path
    //
    CREATE_PATH (
        ch_input.ncbi
    )

    //
    // MODULE: Run ncbigenomedownlaod
    //

    NCBIGENOMEDOWNLOAD ( 
        CREATE_PATH.out.meta,
        CREATE_PATH.out.accession,
        [],
        params.groups
    )
    ch_versions = ch_versions.mix(NCBIGENOMEDOWNLOAD.out.versions.first())
    
    //
    // MODULE: Run GFFREAD
    //

    GFFREAD ( 
        NCBIGENOMEDOWNLOAD.out.fna.mix( ch_input.local.map { [it[0],file(it[1])] } ),
        NCBIGENOMEDOWNLOAD.out.gff.mix( ch_input.local.map { [it[0],file(it[2])] } ) 
    )

    //
    // SUBWORKFLOW: Run TIDK
    //    

    FASTA_EXPLORE_SEARCH_PLOT_TIDK (
        NCBIGENOMEDOWNLOAD.out.fna.mix( ch_input.local.map { [it[0],file(it[1])] } ), []
    )


    //
    // MODULE: Run Orthofinder
    //

    //ortho_ch = [ 
    //    [id:"orthofinder"],
    //    [GFFREAD.out.longest.collect()]
    //]

    ortho_ch = GFFREAD.out.longest.collect().map { it-> [[id:"orthofinder"], it] }

    //ortho_ch.view()

    ORTHOFINDER (
        ortho_ch,
        [[],[]]
    )
    ch_versions = ch_versions.mix(ORTHOFINDER.out.versions)

    //
    // MODULE: Run BUSCO
    //

    BUSCO_BUSCO (  
        GFFREAD.out.proteins_busco, 
        params.busco_mode,
        params.busco_lineage,
        params.busco_lineages_path ?: [],
        params.busco_config ?: []
    )
    ch_versions = ch_versions.mix(BUSCO_BUSCO.out.versions.first())

    //
    // MODULE: Run FastQC
    //
//    FASTQC (
//        ch_samplesheet
//    )
//    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
//    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    //gff            = NCBIGENOMEDOWNLOAD.out.gff
    //fasta          = NCBIGENOMEDOWNLOAD.out.fna
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
