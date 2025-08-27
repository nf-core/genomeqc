/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MERYL_UNIONSUM                       } from '../modules/nf-core/meryl/unionsum/main'
include { MERYL_COUNT                          } from '../modules/nf-core/meryl/count/main'
include { MERQURY_MERQURY                      } from '../modules/nf-core/merqury/merqury/main'
include { CREATE_PATH                          } from '../modules/local/create_path'
include { NCBIGENOMEDOWNLOAD                   } from '../modules/nf-core/ncbigenomedownload/main'
include { PIGZ_UNCOMPRESS as UNCOMPRESS_FASTA  } from '../modules/nf-core/pigz/uncompress/main'
include { PIGZ_UNCOMPRESS as UNCOMPRESS_GXF    } from '../modules/nf-core/pigz/uncompress/main'
include { GENOME_ONLY                          } from '../subworkflows/local/genome_only'
include { GENOME_AND_ANNOTATION                } from '../subworkflows/local/genome_and_annotation'
include { MULTIQC                              } from '../modules/nf-core/multiqc/main'
include { TREE_SUMMARY as TREE_SUMMARY_GENO_ANNO } from '../modules/local/tree_summary'
include { TREE_SUMMARY as TREE_SUMMARY_GENO      } from '../modules/local/tree_summary'
include { paramsSummaryMap                     } from 'plugin/nf-validation'
include { paramsSummaryMultiqc                 } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML               } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText               } from '../subworkflows/local/utils_nfcore_genomeqc_pipeline'
include { validateInputSamplesheet             } from '../subworkflows/local/utils_nfcore_genomeqc_pipeline'
include { FASTA_EXPLORE_SEARCH_PLOT_TIDK       } from '../subworkflows/nf-core/fasta_explore_search_plot_tidk/main'
include { DECONTAMINATION                      } from '../subworkflows/local/decontamination'
include { FCSGX_FETCHDB                        } from '../modules/nf-core/fcsgx/fetchdb/main'
include { BUSCO_SEQS as BUSCO_SEQS_GENOME_ANNO } from '../modules/local/buscos_seqs/main'
include { BUSCO_SEQS as BUSCO_SEQS_GENOME      } from '../modules/local/buscos_seqs/main'
include { SHINY_APP as SHINY_APP_GENOME_ANNO   } from '../modules/local/shiny_app/main'
include { SHINY_APP as SHINY_APP_GENOME        } from '../modules/local/shiny_app/main'

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

    ch_input = ch_samplesheet
                | map {
                    validateInputSamplesheet(it) // Input validation (check local subworkflow for how function works)
                }
                | branch { rows ->
                    ncbi  : rows.size == 3 // channel: [ val(meta), val(refseq), val(fastq) ]
                    local : rows.size == 4 // channel: [ val(meta), val(fasta), val(gxf), val(fastq) ]
                }

    // MODULE: Run create_path

    // ch_input.ncbi is now a 3-element tuple, last element is the fastq.
    // We need to remove it before CREATE_PATH
    ch_input.ncbi
        | map { meta, refseq, fq -> tuple( meta, refseq ) }
        | CREATE_PATH

    // For NCBIGENOMEDOWNLOAD

    ch_ncbi_input = CREATE_PATH.out.accession
                    | multiMap {
                        meta, accession ->
                            meta      : meta
                            accession : accession
                    }

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
    // Perpare fasta channels
    //

    // fasta. We use mix() here becuase when local files are present,
    // then RefSeq IDs should be missing, and viceversa
    fasta        = ch_input.local
                 | map { meta, fasta, gxf, fq -> tuple( meta, fasta) }
                 | mix ( NCBIGENOMEDOWNLOAD.out.fna )

    // Filter fasta files by extension and create channels for each file type
    gz_fasta     = fasta.filter { meta, fasta -> fasta.name.endsWith(".gz") }
    non_gz_fasta = fasta.filter { meta, fasta -> !fasta.name.endsWith(".gz") }

    // Run module uncompress_fasta and combine channels back
    // together so that all the uncompressed files are in channels
    UNCOMPRESS_FASTA ( gz_fasta )
    ch_fasta     = UNCOMPRESS_FASTA.out.file.mix(non_gz_fasta)
    ch_versions  = ch_versions.mix(UNCOMPRESS_FASTA.out.versions.first())

    //
    // Perpare gxf channels
    //

    // gxf. We use mix() here becuase when local files are present,
    // then RefSeq IDs should be missing, and viceversa
    gxf         = ch_input.local
                | map { meta, fasta, gxf, fq ->  tuple( meta,  gxf) }
                | mix ( NCBIGENOMEDOWNLOAD.out.gff )

    // Filter gxf files by extension and create channels for each file type
    gz_gxf      = gxf.filter { meta, gxf -> gxf  && gxf.name.endsWith(".gz")  } // Filter non empty and compressed gxf (channel to be uncompressed)
    non_gz_gxf  = gxf.filter { meta, gxf -> !gxf || !gxf.name.endsWith(".gz") } // Filter empty and uncompressed gxf (not uncompressed)

    // Run module uncompress_GXF and combine channels back
    // together so that all the uncompressed files are in channels
    UNCOMPRESS_GXF( gz_gxf )
    ch_gxf      = UNCOMPRESS_GXF.out.file.mix(non_gz_gxf)
    ch_versions = ch_versions.mix(UNCOMPRESS_GXF.out.versions.first())

    //
    // Perpare gxf channels
    //

    // FASTQ file is optional in the samplesheet.
    // First, get it like you do for gxf and fasta

    ch_fastq = ch_input.ncbi
                | map{ meta, refseq, fq -> tuple( meta, fq ) }
                | mix( ch_input.local.map { meta, fasta, gxf, fq -> tuple( meta, fq ) } )

    //
    // Define multi-channel objects for every process/subworkflow
    //

    // Combine both fasta, gxf and fastq channels into a single multi-channel object
    // using multiMap, so that they are in sync
    // If element (fasta, gxf, fq) is empty, it will return an empty (null) channel
    // Check multimapChannel function below

    ch_input       = ch_fasta // channel: [ val(meta), val(fasta), val(gxf), val(fastq) ]
                   | join(ch_gxf, remainder: true) // reminder: if gxf is null (only necessary for ncbi genomes)
                   | join(ch_fastq, remainder: true)

    // Split into two channels according to the presence/absence of an annotation
    ch_input_anno  = ch_input.filter { meta, fasta, gxf, fastq ->  gxf } // gxf is present. Channel will run on genome and annotation
                   | multimapChannel // Notice only fasta channel and gxf are necessary here
    ch_input_geno  = ch_input.filter { meta, fasta, gxf, fastq ->  !gxf }// gxf is missing. Channel will run on genome only
                   | multimapChannel // Notice only fasta channel is necessary here

    // Merqury
    ch_input_merq  = ch_input.filter { meta, fasta, gxf, fastq -> fastq } // filter rows where fastq is present
                   | multimapChannel // Notice only fasta and fastq channels are necessary here

    // Decontamination subworkflow
    ch_input_decon = ch_fasta.filter { meta, fasta -> meta.taxid } // filter rows where taxid is present. Run decon on those

    // For TIDK the ch_fasta channel will work

    //
    // Run DECONTAMINATION
    //

    // If statement in case people give taxids but no database.
    // This way subworkflow won't try to run (otherwise it'll just fail)
    // Add warning in parameter/input validation plugin
    if ( params.gxdb || params.gxdb_manifiest ) {
        DECONTAMINATION (
            ch_input_decon,
            params.gxdb ?: [],
            params.gxdb_manifiest ?: []
        )
        ch_versions = ch_versions.mix(DECONTAMINATION.out.versions.first())
    }

    //
    // Run TIDK
    //

    if (!params.skip_tidk) {
        FASTA_EXPLORE_SEARCH_PLOT_TIDK (
            ch_fasta,
            []
        )
    }

    // Merqury: Evaluate genome assemblies with k-mers and more
    // https://github.com/marbl/merqury
    // Only run if not skipping and fastq is provided in the samplesheet
        // MODULE: MERYL_COUNT
    MERYL_COUNT(
        ch_input_merq.fq,
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
    ch_versions    = ch_versions.mix(MERYL_UNIONSUM.out.versions.first())
    // MODULE: MERQURY_MERQURY
    ch_merqury_inputs = ch_meryl_union.join(ch_input_merq.fasta)

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

    // Run genome only or genome + gxf
    if (params.genome_only) { // Does this work or should I remove it?
        GENOME_ONLY (
            ch_input.fasta
        )
        ch_multiqc_files = ch_multiqc_files
                         | mix(GENOME_ONLY.out.quast_results.map { meta, results -> results })
                         | mix(GENOME_ONLY.out.busco_short_summaries.map { meta, txt -> txt })
        ch_versions      = ch_versions.mix(GENOME_ONLY.out.versions)
    } else {
        GENOME_ONLY (
            ch_input_geno.fasta
        )
        GENOME_AND_ANNOTATION (
            ch_input_anno.fasta,
            ch_input_anno.gxf
        )
        ch_multiqc_files = ch_multiqc_files
                         | mix(GENOME_AND_ANNOTATION.out.quast_results.map { meta, results -> results })
                         | mix(GENOME_AND_ANNOTATION.out.busco_short_summaries.map { meta, txt -> txt })
        ch_versions      = ch_versions.mix(GENOME_AND_ANNOTATION.out.versions)

        // Number of sequences with more than x complete single copy buscos
        // this should depend on whether protein mode was used or not, not

        BUSCO_SEQS_GENOME_ANNO(
            GENOME_AND_ANNOTATION.out.buscos_per_seqs.map { tables -> [[id:"tables"], tables] }
        )

        BUSCO_SEQS_GENOME(
            GENOME_ONLY.out.buscos_per_seqs.map { tables -> [[id:"tables"], tables] }
        )

        // Prepare channels for tree plot
        ch_tree_genome_anno = GENOME_AND_ANNOTATION.out.tree_data
                            | concat(BUSCO_SEQS_GENOME_ANNO.out.table.map { meta, table -> table})
                            | collect
        ch_tree_genome      = GENOME_ONLY.out.tree_data
                            | concat(BUSCO_SEQS_GENOME.out.table.map { meta, table -> table})
                            | collect

        //
        // MODULE: Run TREE SUMMARY
        //

        TREE_SUMMARY_GENO_ANNO (
            GENOME_AND_ANNOTATION.out.orthofinder,
            ch_tree_genome_anno
        )
        ch_versions      = ch_versions.mix(TREE_SUMMARY_GENO_ANNO.out.versions)

        TREE_SUMMARY_GENO (
            GENOME_ONLY.out.orthofinder,
            ch_tree_genome
        )
        ch_versions      = ch_versions.mix(TREE_SUMMARY_GENO.out.versions)

        //
        // MODULE: Run SHINY APP
        //

        // Prepare script with functions channel
        ch_functions = Channel.fromPath("$projectDir/bin/tree_functions.R", checkIfExists: true)
        ch_app       = Channel.fromPath("$projectDir/bin/shiny_app.R", checkIfExists: true)

        // For genome and annotation
        SHINY_APP_GENOME_ANNO (
            TREE_SUMMARY_GENO_ANNO.out.tables,
            TREE_SUMMARY_GENO_ANNO.out.tree,
            ch_functions,
            ch_app
        )
        ch_versions      = ch_versions.mix(SHINY_APP_GENOME_ANNO.out.versions)

        // For genome only
        SHINY_APP_GENOME (
            TREE_SUMMARY_GENO.out.tables,
            TREE_SUMMARY_GENO.out.tree,
            ch_functions,
            ch_app
        )
        ch_versions      = ch_versions.mix(SHINY_APP_GENOME.out.versions)
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


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def multimapChannel(input) {
   multi_ch = input
            | multiMap {
                meta, fasta, gxf, fq ->
                    fasta : fasta ? tuple( meta, file(fasta) ) : null
                    gxf   : gxf   ? tuple( meta, file(gxf) )   : null
                    fq    : fq    ? tuple( meta, file(fq) )    : null // Only this one is necessary
            }
    return multi_ch
}
