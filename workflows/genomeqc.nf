/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MERYL_UNIONSUM                         } from '../modules/nf-core/meryl/unionsum/main'
include { MERYL_COUNT                         } from '../modules/nf-core/meryl/count/main'
include { MERQURY_MERQURY                     } from '../modules/nf-core/merqury/merqury/main'
include { CREATE_PATH                         } from '../modules/local/create_path'
include { NCBIGENOMEDOWNLOAD                  } from '../modules/nf-core/ncbigenomedownload/main'
include { PIGZ_UNCOMPRESS as UNCOMPRESS_FASTA } from '../modules/nf-core/pigz/uncompress/main'
include { PIGZ_UNCOMPRESS as UNCOMPRESS_GFF   } from '../modules/nf-core/pigz/uncompress/main'
include { GENOME                              } from '../subworkflows/local/genome'
include { GENOME_AND_ANNOTATION               } from '../subworkflows/local/genome_and_annotation'
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
    ch_tree_data = Channel.empty()

    ch_samplesheet
        .map {
            validateInputSamplesheet(it) // Input validation (check local subworkflow)
        }
        .branch {
            ncbi: it.size() == 3
            local: it.size() == 4
        }
        .set { ch_input }

    // MODULE: Run create_path
    // ch_input.ncbi is now a 3-element tuple, last element is the fastq. We need to remove it before CREATE_PATH
    ch_input.ncbi
        | map { [it[0], it[1]] }
        | CREATE_PATH


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
    gff   = NCBIGENOMEDOWNLOAD.out.gff.mix( ch_input.local.map { [it[0],file(it[1])] } )

    // Uncompress files if necessary | Consider using brances as an alternative
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

    // FASTQ file is optional in the samplesheet. 
    // First, get it like you do for gff and fasta
    ch_fastq = ch_input.local.map{ [it[0], it[2]] }.mix(ch_input.ncbi.map{ [it[0],it[2]] })
    // Then, check to see that element 1 is not empty, and if not, make it file()
    // You have to do this because if you pass in file() in the initial map, 
    // it'll fail if you don't supply a fastq, because you can't pass an empty to file()
    ch_fastq
        | map{meta, fq -> fq ? [meta, file(fq)] : [meta, fq]}
        | filter { meta, fq -> fq && fq.name =~ /(\.fastq|\.fq|\.fastq\.gz|\.fq\.gz)$/ }
        | set {ch_fastq}
    
    //
    // Run TIDK
    //
    FASTA_EXPLORE_SEARCH_PLOT_TIDK (
        ch_fasta,
        []
    )

    // Merqury: Evaluate genome assemblies with k-mers and more
    // https://github.com/marbl/merqury
    // Only run if not skipping and fastq is provided in the samplesheet
    if (!params.merqury_skip && ch_fastq) {
        // MODULE: MERYL_COUNT
        MERYL_COUNT(
            ch_fastq,
            params.kvalue 
        )
        ch_meryl_db = MERYL_COUNT.out.meryl_db
        ch_versions = ch_versions.mix(MERYL_COUNT.out.versions.first())
        // MODULE: MERYL_UNIONSUM
        MERYL_UNIONSUM(
            ch_meryl_db,
            params.kvalue
        )
        ch_meryl_union = MERYL_UNIONSUM.out.meryl_db
        ch_versions = ch_versions.mix(MERYL_UNIONSUM.out.versions.first())
        // MODULE: MERQURY_MERQURY
        ch_meryl_union
            | join(ch_fastq)
            | set {ch_merqury_inputs}
        MERQURY_MERQURY ( ch_merqury_inputs )
        ch_merqury_qv                           = MERQURY_MERQURY.out.assembly_qv
        ch_merqury_stats                        = MERQURY_MERQURY.out.stats
        ch_merqury_spectra_cn_fl_png            = MERQURY_MERQURY.out.spectra_cn_fl_png
        ch_merqury_spectra_asm_fl_png           = MERQURY_MERQURY.out.spectra_asm_fl_png
        ch_hapmers_blob_png                     = MERQURY_MERQURY.out.hapmers_blob_png
        ch_merqury_outputs                      = ch_merqury_qv
                                                | mix(ch_merqury_stats)
                                                | mix(ch_merqury_spectra_cn_fl_png)
                                                | mix(ch_merqury_spectra_asm_fl_png)
                                                | mix(ch_hapmers_blob_png)
                                                | flatMap { meta, data -> data }
        ch_versions                             = ch_versions.mix(MERQURY_MERQURY.out.versions.first())
    }


    // Run genome only or genome + gff
    if (params.genome_only) {
        GENOME (
            ch_fasta
        )
        ch_multiqc_files = ch_multiqc_files.mix(GENOME.out.quast.collect())
        ch_multiqc_files = ch_multiqc_files.mix(GENOME.out.busco.collect())
    } else {
        GENOME_AND_ANNOTATION (
            ch_fasta,
            ch_gff
        )
        ch_multiqc_files = ch_multiqc_files.mix(GENOME_AND_ANNOTATION.out.quast.collect())
        ch_multiqc_files = ch_multiqc_files.mix(GENOME_AND_ANNOTATION.out.busco.collect())
    }

    ch_multiqc_files.view()

    //
    // MODULE: Run TREE SUMMARY
    //  
    TREE_SUMMARY (
        GENOME_AND_ANNOTATION.out.orthofinder,
        GENOME_AND_ANNOTATION.out.tree_data
    )


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
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
