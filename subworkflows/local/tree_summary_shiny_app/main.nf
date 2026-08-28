include { TREESUMMARY as TREE_SUMMARY_GENO_ANNO    } from '../../../modules/local/treesummary/main'
include { TREESUMMARY as TREE_SUMMARY_GENO         } from '../../../modules/local/treesummary/main'
include { SHINYAPP as SHINY_APP_GENOME_ANNO        } from '../../../modules/local/shinyapp/main'
include { SHINYAPP as SHINY_APP_GENOME             } from '../../../modules/local/shinyapp/main'


workflow TREE_SUMMARY_SHINY_APP {

    take:
    ch_geno_anno_busco_summary
    ch_prot_anno_busco_summary
    ch_geno_busco_summary
    ch_tree_inputs_genome_anno
    ch_tree_inputs_genome
    ch_geno_anno_orthofinder
    ch_geno_orthofinder
    ch_te_table
    ch_fcs_table
    val_skip_busco       // val: boolean - skip genome-only tree summary/shiny app (needs BUSCO markers)
    val_container_engine // val: container engine used to launch the interactive Shiny app - 'docker' or 'podman'

    main:

    //
    // MODULE: Run TREE SUMMARY
    //

    // Prepare busco channel for genome and annotation
    // First for genome completeness
    ch_busco_geno_anno1 = ch_geno_anno_busco_summary
                        | map { _meta, file -> file }
                        | collect
                        | map { files -> tuple( [id:"busco_geno_anno"], files )}
                        | ifEmpty([[],[]])

    // Then for annotation completeness
    ch_busco_geno_anno2 = ch_prot_anno_busco_summary
                        | map { _meta, file -> file }
                        | collect
                        | map { files -> tuple( [id:"busco_geno_anno"], files )}
                        | ifEmpty([[],[]])

    // Combine both channels into a multi-channel object
    ch_busco_geno_anno  = ch_busco_geno_anno1.join(ch_busco_geno_anno2)
                        | multiMap {
                            meta, geno_files, prot_files ->
                                geno      : geno_files ? tuple( meta, geno_files ) : [[],[]]
                                prot      : prot_files ? tuple( meta, prot_files ) : [[],[]]
                        }

    // Prepare busco channel for genome only
    ch_busco_geno       = ch_geno_busco_summary
                        | map { _meta, file -> file }
                        | collect
                        | map { files -> tuple( [id:"busco_geno"], files )}
                        | ifEmpty([[],[]])

    // Prepare channels for shiny app
    ch_functions = channel.fromPath("$projectDir/bin/tree_functions.R", checkIfExists: true)
    ch_app       = channel.fromPath("$projectDir/bin/shiny_app.R", checkIfExists: true)

    // Run TREE SUMMARY for genome and annotation
    TREE_SUMMARY_GENO_ANNO (
        ch_geno_anno_orthofinder,
        ch_busco_geno_anno.geno,
        ch_busco_geno_anno.prot,
        ch_te_table,
        ch_fcs_table,
        ch_tree_inputs_genome_anno
    )

    // Run SHINY APP for genome and annotation
    SHINY_APP_GENOME_ANNO (
        TREE_SUMMARY_GENO_ANNO.out.tables.join(TREE_SUMMARY_GENO_ANNO.out.tree, by:0),
        ch_functions,
        ch_app,
        val_container_engine
    )

    if (!val_skip_busco) { // Tree sumary on genome only is only run if BUSCO is not skipped (markers needed for tree summary)
        // Run TREE SUMMARY for genome only
        TREE_SUMMARY_GENO (
            ch_geno_orthofinder,
            ch_busco_geno,
            [[],[]], // No busco proteins for genome only (busco runs on genome)
            ch_te_table,
            ch_fcs_table,
            ch_tree_inputs_genome
        )
        // Run SHINY APP for genome only
        SHINY_APP_GENOME (
            TREE_SUMMARY_GENO.out.tables.join(TREE_SUMMARY_GENO.out.tree, by:0),
            ch_functions,
            ch_app,
            val_container_engine
        )

    }

    emit:
    geno_anno_tables = TREE_SUMMARY_GENO_ANNO.out.tables
    geno_anno_tree   = TREE_SUMMARY_GENO_ANNO.out.tree
    geno_tables      = !val_skip_busco ? TREE_SUMMARY_GENO.out.tables      : channel.empty()
    geno_tree        = !val_skip_busco ? TREE_SUMMARY_GENO.out.tree        : channel.empty()
}
