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
        | map {
            validateInputSamplesheet(it) // Input validation (check local subworkflow)
        }
        | branch {
            ncbi  : it.size() == 3
            local : it.size() == 4
        }
        | set { ch_input }

    // MODULE: Run create_path

    // ch_input.ncbi is now a 3-element tuple, last element is the fastq.
    // We need to remove it before CREATE_PATH
    ch_input.ncbi
        | map { meta, refseq, fq -> tuple( meta, refseq ) }
        | CREATE_PATH

    // For NCBIGENOMEDOWNLOAD

    CREATE_PATH.out.accession
        | multiMap {
            meta, accession ->
                meta      : meta
                accession : accession
        }
        | set { ch_ncbi_input }

    //
    // MODULE: Run ncbigenomedownlaod for RefSeq IDs
    //

    NCBIGENOMEDOWNLOAD ( 
        ch_ncbi_input.meta,
        ch_ncbi_input.accession,
        [],
        params.groups
    )
    ch_versions = ch_versions.mix(NCBIGENOMEDOWNLOAD.out.versions.first())
    
    //
    // Define gff and fasta channels
    //

    //fasta = NCBIGENOMEDOWNLOAD.out.fna.mix( ch_input.local.map { meta, fasta, gff, fq -> tuple( meta, file(fasta) ) } )
    //gff   = NCBIGENOMEDOWNLOAD.out.gff.mix( ch_input.local.map { meta, fasta, gff, fq -> tuple( meta, file(gff) ) } )
    
    // gff. We use mix() here becuase when local files are present,
    // then RefSeq IDs should be missing, and viceversa
    ch_input.local
        | map { meta, fasta, gff, fq -> tuple( meta, file(fasta) ) }
        | mix ( NCBIGENOMEDOWNLOAD.out.fna )
        | set { fasta }
    // fasta. We use mix() here becuase when local files are present, then RefSeq IDs should be missing, and viceversa
    ch_input.local
        | map { meta, fasta, gff, fq -> tuple( meta, file(gff) ) }
        | mix ( NCBIGENOMEDOWNLOAD.out.gff )
        | set { gff }


    // Filter fasta files by extension and create channels for each file type
    gz_fasta = fasta.filter { it[1].name.endsWith(".gz") }
    non_gz_fasta = fasta.filter { !it[1].name.endsWith(".gz") }

    // Run module uncompress_fasta

    UNCOMPRESS_FASTA (gz_fasta )
    ch_versions = ch_versions.mix(UNCOMPRESS_FASTA.out.versions.first())

    // Filter gff files by extension and create channels for each file type

    gz_gff = gff.filter { it[1].name.endsWith(".gz") }
    non_gz_gff = gff.filter { !it[1].name.endsWith(".gz") }

    // Run module uncompress_GFF

    UNCOMPRESS_GFF(gz_gff)
    ch_versions = ch_versions.mix(UNCOMPRESS_GFF.out.versions.first())

    // Combine the channels back together so that all the uncompressed files are in channels 

    ch_fasta  = UNCOMPRESS_FASTA.out.file.mix(non_gz_fasta)
    ch_gff = UNCOMPRESS_GFF.out.file.mix(non_gz_gff)

    //
    // Define fastq input channel
    //

    // FASTQ file is optional in the samplesheet. 
    // First, get it like you do for gff and fasta

    ch_input.ncbi
        | map{ meta, refseq, fq -> tuple( meta, fq ) }
        | mix( ch_input.local.map{ meta, fasta, gff, fq -> tuple( meta, fq ) } )
        | set { ch_fastq }

    // Then, check to see that element 1 is not empty, and if not, make it file()
    // You have to do this because if you pass in file() in the initial map, 
    // it'll fail if you don't supply a fastq, because you can't pass an empty to file()

    ch_fastq
        | map{meta, fq -> fq ? [meta, file(fq)] : [meta, fq]}
        | filter { meta, fq -> fq && fq.name =~ /(\.fastq|\.fq|\.fastq\.gz|\.fq\.gz)$/ }
        | set {ch_fastq}
    
    //
    // Define multi-channel object
    //

    // Combine both fasta, gff and fastq channels into a single multi-channel object using multiMap, so that they are in sync all the time

    ch_fasta
        | combine(ch_gff, by:0) // by:0 | Only combine when both channels share the same id
        | combine(ch_fastq, by:0)
        | multiMap {
            meta, fasta, gff, fastq ->
                fasta : tuple( meta, fasta )
                gff   : tuple( meta, gff )
                fq    : tuple( meta, fastq )
        }
        | set { ch_input }

    //
    // Run TIDK
    //
    
    if (!params.skip_tidk) {
        FASTA_EXPLORE_SEARCH_PLOT_TIDK (
            ch_input.fasta,
            []
        )
    }

    // Merqury: Evaluate genome assemblies with k-mers and more
    // https://github.com/marbl/merqury
    // Only run if not skipping and fastq is provided in the samplesheet
    if (!params.merqury_skip && ch_input.fq) {
        // MODULE: MERYL_COUNT
        MERYL_COUNT(
            ch_input.fq,
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
            | join(ch_input.fasta)
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
            ch_input.fasta
        )
    } else {
        GENOME_AND_ANNOTATION (
            ch_input.fasta,
            ch_input.gff
        )
        
        //
        // MODULE: Run TREE SUMMARY
        //  

        TREE_SUMMARY (
            GENOME_AND_ANNOTATION.out.orthofinder,
            GENOME_AND_ANNOTATION.out.tree_data
        )
    }

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
