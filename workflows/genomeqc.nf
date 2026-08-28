/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTA_QV_MERQURY                         } from '../subworkflows/local/fasta_qv_merqury/main'
include { INPUT_PREPARATION                        } from '../subworkflows/local/input_preparation/main'
include { GENOME_ONLY                              } from '../subworkflows/local/genome_only/main'
include { GENOME_AND_ANNOTATION                    } from '../subworkflows/local/genome_and_annotation/main'
include { TREE_SUMMARY_SHINY_APP                   } from '../subworkflows/local/tree_summary_shiny_app/main'
include { FASTA_EXPLORE_SEARCH_PLOT_TIDK           } from '../subworkflows/nf-core/fasta_explore_search_plot_tidk/main'
include { DECONTAMINATION                          } from '../subworkflows/local/decontamination/main'
include { TABLE_FCSGX                              } from '../modules/local/table/fcsgx/main'
include { BUSCO_SEQS as BUSCO_SEQS_GENOME_ANNO     } from '../modules/local/busco/seqs/main'
include { BUSCO_SEQS as BUSCO_SEQS_GENOME          } from '../modules/local/busco/seqs/main'
include { REPORT_HTML                              } from '../modules/local/report/html/main'
include { REPORT_EXCEL                             } from '../modules/local/report/excel/main'
include { FASTA_ANNOTATE_TE                        } from '../subworkflows/local/fasta_annotate_te/main'
include { MULTIQC                                  } from '../modules/nf-core/multiqc/main'
include { validateInputSamplesheet                 } from '../subworkflows/local/utils_nfcore_genomeqc_pipeline'
include { multimapChannel                          } from '../subworkflows/local/utils_nfcore_genomeqc_pipeline'
include { paramsSummaryMap                         } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                     } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                   } from '../subworkflows/local/utils_nfcore_genomeqc_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GENOMEQC {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir
    val_groups               // val: NCBI genome download assembly groups (e.g. 'all')
    val_busco_lineages_path  // val: path to a pre-staged BUSCO lineages directory, or null
    val_busco_lineage        // val: BUSCO lineage name (e.g. 'hymenoptera_odb10'), 'auto', or null
    val_gxdb                 // val: path to a local FCS-GX database, or []
    val_gxdb_manifest        // val: path to an FCS-GX database manifest, or []
    val_ramdisk              // val: path to a ramdisk for FCS-GX, or []
    val_repeat               // val: telomeric repeat motif for tidk, or null
    val_skip_tidk            // val: boolean - skip the tidk subworkflow
    val_kvalue               // val: Merqury k-mer size
    val_RM_download_db       // val: boolean - download a DFAM RepeatMasker database partition
    val_RM_db                // val: list of DFAM database partition names to download
    val_famdb_library        // val: path (or glob) to a pre-staged famdb library, or null
    val_famdb_lineage        // val: lineage string for famdb extraction (e.g. 'hymenoptera'), or null
    val_run_repeatmodeler    // val: boolean - run de novo RepeatModeler (slow, adds 24-48 h per genome)
    val_te_clusterer         // val: clustering tool - 'linclust' (default), 'mmseqs', or 'cdhit'
    val_te                   // val: 'hite' or 'repeatmasker', or null to skip TE annotation
    val_skip_busco           // val: boolean - skip BUSCO (and everything downstream of it)
    val_busco_config         // val: path to a BUSCO config file, or []
    val_busco_clean          // val: boolean - clean up intermediate BUSCO files, or []
    val_ortho_version        // val: OrthoFinder version - 'v2' or 'v3'
    val_validation_tool      // val: GXF validation/standardisation tool - 'agat' or 'gffread'
    val_container_engine     // val: container engine used to launch the interactive Shiny app - 'docker' or 'podman'

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()

    ch_input = ch_samplesheet
                | map { samplesheet ->
                    validateInputSamplesheet(samplesheet) // Input validation (check local subworkflow for how function works)
                }
                | branch { rows ->
                    ncbi  : rows.size == 3 // channel: [ val(meta), val(refseq), val(fastq) ]
                    local : rows.size == 4 // channel: [ val(meta), val(fasta), val(gxf), val(fastq) ]
                }

    //
    // SUBWORKFLOW: Prepare input files (fasta, gxf, fastq) and BUSCO database
    //

    INPUT_PREPARATION (
        ch_input.ncbi,
        ch_input.local,
        val_groups,
        val_busco_lineages_path,
        val_busco_lineage
    )

    // ch_fasta contains ALL samples regardless of annotation presence
    ch_fasta       = INPUT_PREPARATION.out.fasta
    ch_input_decon = INPUT_PREPARATION.out.input_decon
    ch_busco_db    = INPUT_PREPARATION.out.busco_db

    // Prepare channels for genome only, genome + annotation, and Merqury subworkflows
    ch_input_anno  = INPUT_PREPARATION.out.input
                   | filter { _meta, _fasta, gxf, _fastq ->  gxf }
                   | multimapChannel
    // Not redundant, not the same as ch_fasta
    ch_input_geno  = INPUT_PREPARATION.out.input
                   | filter { _meta, _fasta, gxf, _fastq -> !gxf }
                   | multimapChannel
    ch_input_merq  = INPUT_PREPARATION.out.input
                   | filter { _meta, _fasta, _gxf, fastq -> fastq }
                   | multimapChannel

    //
    // SUBWORKFLOW: Run DECONTAMINATION
    //

    // If statement in case people give taxids but no database.
    // This way subworkflow won't try to run (otherwise it'll just fail)
    // Add warning in parameter/input validation plugin
    ch_fcs_table = channel.empty()
    if ( val_gxdb || val_gxdb_manifest ) {
        DECONTAMINATION (
            ch_input_decon,
            val_ramdisk ?: [],
            val_gxdb ?: [],
            val_gxdb_manifest ? file(val_gxdb_manifest) : []
        )

        // Parse each species' FCS-GX report into a combined contamination-summary
        // table for the tree plot (this is only needed for the tree plot, which
        // requires a single TSV table, same as the TE-composition table above)
        TABLE_FCSGX (
            DECONTAMINATION.out.fcs_gx_report
                .map { _meta, f -> f }
                .collect()
                .map { f -> tuple([id:'fcs_table'], f) }
        )
        ch_fcs_table = TABLE_FCSGX.out.table
    }

    //
    // SUBWORKFLOW: Run TIDK
    //

    ch_repeat = val_repeat ? ch_fasta.map { meta, _fasta_repeat -> [ meta, val_repeat ] } : channel.empty()

    if (!val_skip_tidk) {
        FASTA_EXPLORE_SEARCH_PLOT_TIDK (
            ch_fasta,
            ch_repeat
        )
    }

    //
    // SUBWORKFLOW: Run Merqury
    //

    FASTA_QV_MERQURY (
        ch_input_merq.fasta,
        ch_input_merq.fq,
        val_kvalue
    )

    //
    // SUBWORKFLOW: Annotate transposable elements
    //

    ch_te_table = channel.empty()
    ch_te_table_collect = channel.empty()

    // Skip download when a pre-staged famdb library is already provided
    ch_rm_db_input = (val_RM_download_db && val_RM_db && !val_famdb_library)
                   ? channel.fromList(val_RM_db)
                       .map { db -> tuple([id: file(db).getBaseName()], db) }
                   : channel.empty()

    // Accepts a single file path or a glob pattern (e.g. '/path/FamDB*')
    ch_famdb_lib_input = val_famdb_library
                       ? channel.fromPath(val_famdb_library)
                           .map { path -> tuple([id: path.baseName], path) }
                       : channel.empty()

    if (val_te) {
        FASTA_ANNOTATE_TE (
            ch_fasta,
            ch_rm_db_input.ifEmpty([[],[]]),
            ch_famdb_lib_input.ifEmpty([[],[]]),
            val_famdb_lineage ?: '',
            val_run_repeatmodeler,
            val_te_clusterer,
            val_te
        )
        ch_te_table = FASTA_ANNOTATE_TE.out.tbl_tsv
        ch_te_table_collect = FASTA_ANNOTATE_TE.out.tbl_collected
    }

    //
    // SUBWORKFLOW: Run genome only or genome + annotation subworkflows to get QC metrics
    //

    GENOME_ONLY (
        ch_input_geno.fasta,
        ch_busco_db.ifEmpty([]),
        val_skip_busco,
        val_busco_lineage,
        val_busco_config,
        val_busco_clean,
        val_ortho_version
    )
    ch_geno_busco_summary = GENOME_ONLY.out.busco_short_summaries
    ch_geno_orthofinder   = GENOME_ONLY.out.orthofinder

    GENOME_AND_ANNOTATION (
        ch_input_anno.fasta,
        ch_input_anno.gxf,
        ch_busco_db.ifEmpty([]),
        val_validation_tool,
        val_ortho_version,
        val_skip_busco,
        val_busco_lineage,
        val_busco_config,
        val_busco_clean
    )

    // Define channels for tree summary and shiny app subworkflow
    ch_geno_anno_busco_summary = GENOME_AND_ANNOTATION.out.busco_short_summaries_geno
    ch_prot_anno_busco_summary = GENOME_AND_ANNOTATION.out.busco_short_summaries_prot
    ch_geno_anno_orthofinder   = GENOME_AND_ANNOTATION.out.orthofinder


    ch_multiqc_files = ch_multiqc_files
                     | mix(GENOME_AND_ANNOTATION.out.quast_results.map { _meta, results -> results })
                     | mix(ch_prot_anno_busco_summary.map { _meta, txt -> txt })

    //
    // MODULE: run BUSCO SEQS
    //

    // Number of sequences with more than x complete single copy buscos
    // this should depend on whether protein mode was used or not
    if(!val_skip_busco){
        BUSCO_SEQS_GENOME_ANNO(
            GENOME_AND_ANNOTATION.out.buscos_per_seqs.map { tables -> [[id:"tables"], tables] }
        )

        BUSCO_SEQS_GENOME(
            GENOME_ONLY.out.buscos_per_seqs.map { tables -> [[id:"tables"], tables] }
        )

        // Prepare channels for tree plot
        ch_tree_genome_anno = GENOME_AND_ANNOTATION.out.tree_data
                            | concat(BUSCO_SEQS_GENOME_ANNO.out.table.map { _meta, table -> table})
                            | collect
        ch_tree_genome      = GENOME_ONLY.out.tree_data
                            | concat(BUSCO_SEQS_GENOME.out.table.map { _meta, table -> table})
                            | collect

        // Optional channel for HTML report: empty list if no busco seqs tables were produced.
        // collectFile merges both modes' tables (disjoint species) into one TSV with a single
        // header, since generate_report.py/generate_excel.py expect exactly one file.
        ch_busco_seqs_ga_file = BUSCO_SEQS_GENOME_ANNO.out.table
                              | mix (BUSCO_SEQS_GENOME.out.table)
                              | map { _meta, f -> f }
                              | collectFile(name: 'n_seqs_above_x_buscos.tsv', keepHeader: true)
                              | ifEmpty([])

    } else { // If BUSCO is not run
        ch_tree_genome_anno   = GENOME_AND_ANNOTATION.out.tree_data.collect()
        ch_tree_genome        = GENOME_ONLY.out.tree_data.collect()
        ch_busco_seqs_ga_file = channel.empty()
    }

    //
    // SUBWORKFLOW: Run static and interactive tree summary plots
    //

    TREE_SUMMARY_SHINY_APP (
        ch_geno_anno_busco_summary,
        ch_prot_anno_busco_summary,
        ch_geno_busco_summary,
        ch_tree_genome_anno,
        ch_tree_genome,
        ch_geno_anno_orthofinder,
        ch_geno_orthofinder,
        ch_te_table.ifEmpty([[],[]]),
        ch_fcs_table.ifEmpty([[],[]]),
        val_skip_busco,
        val_container_engine,
    )

    //
    // MODULE: Run HTML REPORT (genome + annotation mode)
    //
    ch_report_busco = ch_geno_anno_busco_summary
                    | map { _meta, f -> f }
                    | mix( GENOME_ONLY.out.busco_batch_summaries.map { _meta, f -> f } )
                    | collect
                    | first
    // Protein BUSCO only exists for annotated genomes; ifEmpty keeps the
    // report running when no annotations are present in this branch.
    ch_report_busco_prot = ch_prot_anno_busco_summary
                    | map { _meta, f -> f }
                    | collect
                    | ifEmpty( [] )
                    | first
    ch_tidk                = val_skip_tidk ? channel.value([]) : FASTA_EXPLORE_SEARCH_PLOT_TIDK.out.aposteriori_tsv.map { _meta, f -> f }.collect().first()
    ch_report_tidk_apriori = (val_skip_tidk || !val_repeat) ? channel.value([]) : FASTA_EXPLORE_SEARCH_PLOT_TIDK.out.apriori_tsv.map { _meta, f -> f }.collect()

    ch_fcsgx  = (val_gxdb || val_gxdb_manifest) ? DECONTAMINATION.out.fcs_gx_report.map { _meta, f -> f }.collect().first()  : channel.value([])
    ch_fcsadp = (val_gxdb || val_gxdb_manifest) ? DECONTAMINATION.out.adaptor_report.map { _meta, f -> f }.collect().first()  : channel.value([])
    ch_tiara  = (val_gxdb || val_gxdb_manifest) ? DECONTAMINATION.out.tiara_cleaned.map  { _meta, f -> f }.collect().first()  : channel.value([])
    ch_excel_agat = GENOME_AND_ANNOTATION.out.agat_stats.map { _meta, f -> f }.collect().ifEmpty([])
    ch_excel_quast  = GENOME_AND_ANNOTATION.out.quast_tsv.map { _meta, f -> f }
                    | mix( GENOME_ONLY.out.quast_tsv.map { _meta, f -> f } )
                    | collect

    REPORT_HTML (
        ch_report_busco,
        ch_report_busco_prot,
        ch_excel_agat,
        ch_excel_quast.ifEmpty([]),
        ch_tidk,
        ch_report_tidk_apriori,
        ch_fcsgx,
        ch_fcsadp,
        ch_tiara,
        ch_busco_seqs_ga_file.ifEmpty([]),
        ch_te_table_collect.ifEmpty([])
    )

    //
    // MODULE: Run EXCEL REPORT (genome + annotation mode)
    //
    REPORT_EXCEL (
        ch_report_busco,
        ch_report_busco_prot,
        ch_excel_quast,
        ch_excel_agat.ifEmpty([]),
        ch_tidk,
        ch_fcsgx,
        ch_fcsadp,
        ch_tiara,
        ch_busco_seqs_ga_file.ifEmpty([]),
        ch_te_table_collect.ifEmpty([])
    )

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_'  +  'genomeqc_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))

    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'genomeqc'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )

    emit:
    multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
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
