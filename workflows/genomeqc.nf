/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CREATE_PATH                         } from '../modules/local/create_path'
include { NCBIGENOMEDOWNLOAD                  } from '../modules/nf-core/ncbigenomedownload/main'
include { LONGEST                             } from '../modules/local/longest'
include { PIGZ_UNCOMPRESS as UNCOMPRESS_FASTA } from '../modules/nf-core/pigz/uncompress/main'
include { PIGZ_UNCOMPRESS as UNCOMPRESS_GFF   } from '../modules/nf-core/pigz/uncompress/main'
include { BUSCO_BUSCO                         } from '../modules/nf-core/busco/busco/main'
include { GFFREAD                             } from '../modules/nf-core/gffread/main'
include { ORTHOFINDER                         } from '../modules/nf-core/orthofinder/main'
include { FASTQC                              } from '../modules/nf-core/fastqc/main'
include { MULTIQC                             } from '../modules/nf-core/multiqc/main'
include { TREE_SUMMARY                        } from '../modules/local/tree_summary'
include { paramsSummaryMap                    } from 'plugin/nf-validation'
include { paramsSummaryMultiqc                } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML              } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText              } from '../subworkflows/local/utils_nfcore_genomeqc_pipeline'
include { validateInputSamplesheet            } from '../subworkflows/local/utils_nfcore_genomeqc_pipeline'
include { FASTA_EXPLORE_SEARCH_PLOT_TIDK      } from '../subworkflows/nf-core/fasta_explore_search_plot_tidk/main'

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
    // Define gff and fasta varliables
    //

    fasta = NCBIGENOMEDOWNLOAD.out.fna.mix( ch_input.local.map { [it[0],file(it[2])] } )
    gff =   NCBIGENOMEDOWNLOAD.out.gff.mix( ch_input.local.map { [it[0],file(it[1])] } )
    
    // Uncompress fasta if necessary | Consider using brances as an alternative

    if (fasta.map { it[1].endsWith(".gz") } ) {
        ch_fasta = UNCOMPRESS_FASTA ( fasta ).file
    } else {
        ch_fasta = fasta
    }
    
    // Uncompress gff if necessary

    if (gff.map { it[1].endsWith(".gz") } ) {
        ch_gff = UNCOMPRESS_GFF ( gff ).file
    } else {
        ch_gff = gff
    }

    //
    // Run TIDK
    //

    FASTA_EXPLORE_SEARCH_PLOT_TIDK (
        ch_fasta,
        []
    )

    //
    // Run AGAT longest isoform
    //

    LONGEST (
        ch_gff
    )

    //
    // Run uncompressed, GFFREAD requires uncompressed fasta as input
    //
    // Running pigz uncompress should be optional
    
    
    //PIGZ_UNCOMPRESS (
    //    NCBIGENOMEDOWNLOAD.out.fna.map { [it[0],file(it[1])] }.mix( ch_input.local.map { [it[0],file(it[1])] } )
    //)
    
    //
    // Run GFFREAD
    //

    //ch_fasta.view()
    //LONGEST.out.longest_proteins.view()

    ch_long_gff = LONGEST.out.longest_proteins
    

    // Step 1: View the content of each channel before joining
    //ch_long_gff.view()
    //ch_fasta.view()

    inputChannel = ch_long_gff.combine(ch_fasta, by: 0)

    // Split the input channel into two channels
    gffChannel = inputChannel.map { tuple ->
        // Extracting the GFF path and ID
        [tuple[0], tuple[1]]
    }

    fnaChannel = inputChannel.map { tuple ->
        // Extracting only the FNA path
        tuple[2]
    }

    // Usage example: Print the outputs of both channels
    gffChannel.view()
    fnaChannel.view()

    GFFREAD ( 
        gffChannel , fnaChannel
    )

    //
    // MODULE: Run Orthofinder
    //

    ortho_ch = GFFREAD.out.gffread_fasta.map { it[1] }.collect().map { it-> [[id:"orthofinder"], it] }

    ORTHOFINDER (
        ortho_ch,
        [[],[]]
    )
    ch_versions = ch_versions.mix(ORTHOFINDER.out.versions)

    //
    // MODULE: Run BUSCO
    //

    BUSCO_BUSCO (  
        GFFREAD.out.gffread_fasta, 
        params.busco_mode,
        params.busco_lineage,
        params.busco_lineages_path ?: [],
        params.busco_config ?: []
    )
    ch_versions = ch_versions.mix(BUSCO_BUSCO.out.versions.first())

    //
    // MODULE: Run TREE SUMMARY
    //  

    ORTHOFINDER.out.orthofinder.view()

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
