include { QUAST                               } from '../../../modules/nf-core/quast/main'
include { BUSCO_BUSCO as BUSCO_GENOME         } from '../../../modules/nf-core/busco/busco/main'
include { IDEOGRAM_GENOME                     } from '../../../modules/local/ideogram/genome/main'
include { ORTHOFINDER as ORTHOFINDER_V3       } from '../../../modules/nf-core/orthofinder/main'
include { ORTHOFINDERV2                       } from '../../../modules/local/orthofinderv2/main'
include { BUSCO_TSVTOGFF                      } from '../../../modules/local/busco/tsvtogff/main'
include { ORTHOSEQCOUNT                       } from '../../../modules/local/orthoseqcount'
include { GAWK                                } from '../../../modules/nf-core/gawk/main'

workflow GENOME_ONLY {

    take:
    ch_fasta          // channel: [ val(meta), [ fasta ] ]
    ch_busco_db       // channel: path | []
    val_skip_busco    // val: boolean - skip BUSCO (and everything downstream of it: ideogram, orthofinder, orthoseqcount)
    val_busco_lineage // val: BUSCO lineage name (e.g. 'hymenoptera_odb10') or 'auto'
    val_busco_config  // val: path to a BUSCO config file, or []
    val_busco_clean   // val: boolean - clean up intermediate BUSCO files, or []
    val_ortho_version // val: OrthoFinder version - 'v2' or 'v3'

    main:
    ch_fasta.view { meta, _fasta -> "Running ${meta.id} on genome only mode"}


    // For tree plot
    ch_tree_data = channel.empty()

    //
    // MODULE: Run Quast
    //

    QUAST (
        ch_fasta,
        [[],[]],
        [[],[]]
    )
    ch_tree_data = ch_tree_data.mix(QUAST.out.tsv.map { tuple -> tuple[1] })

    if (!val_skip_busco) {
        BUSCO_GENOME (
            ch_fasta,
            "genome", // hardcoded, other options ('proteins', 'transcriptome') make no sense
            val_busco_lineage,
            ch_busco_db,
            val_busco_config ?: [],
            val_busco_clean ?: []
        )
        //ch_tree_data  = ch_tree_data.mix(BUSCO_GENOME.out.batch_summary.collect { meta, file -> file })

        //
        // GAWK
        //
        // Use GAWK to change ID from file name to meta.id

        GAWK (
            BUSCO_GENOME.out.batch_summary,
            [],
            false
        )
        ch_tree_data  = ch_tree_data.mix(GAWK.out.output.collect { _meta, file -> file })

        //
        // BUSCO Ideogram
        //

        ch_full_table = BUSCO_GENOME.out.full_table

        // Combined ch_fasta and BUSCO output channel into a single channel for ideogram
        ch_input_ideo = ch_fasta
                      | combine(ch_full_table, by:0)


        IDEOGRAM_GENOME (
            ch_input_ideo
        )

        //
        // Orthofinder
        //

        // Prepare data
        ch_busco_proteins = BUSCO_GENOME.out.single_copy_faa
                          | flatMap { meta, faas ->
                                         faas.collect { faa -> [meta, file(faa)]  }
                          }
                          | collectFile { meta, faas ->
                                            [ "${meta.id}.fasta", faas ]
                          }
                          | collect
                          | filter { file_paths ->
                                        file_paths.size() >= 3 // Ensure there are enough BUSCO proteins, otherwise skip orthofinder
                          }
                          | map { file_paths ->
                                    [[ id: 'orthofinder', mode: 'genome' ], file_paths.toSorted { a, b -> a.name <=> b.name }] // for deterministic behaviour
                          }

        //Run orthofinder
        if (val_ortho_version == 'v3') {
            ORTHOFINDER_V3 (
                ch_busco_proteins,
                [[],[]]
            )
            ch_orthofinder = ORTHOFINDER_V3.out.orthofinder
        } else if (val_ortho_version == 'v2' ) {
            ORTHOFINDERV2 (
                ch_busco_proteins
            )
            ch_orthofinder = ORTHOFINDERV2.out.orthofinder
        }
        // Transform tsv to gff for orthoseqcount module
        BUSCO_TSVTOGFF (
            BUSCO_GENOME.out.busco_dir
        )

        //
        // MODULE: Run ORTHOSEQCOUNT
        //

        ORTHOSEQCOUNT (
            ch_orthofinder.map { _meta, folder ->
                file("${folder}/Orthogroups/Orthogroups.tsv")
            },
            BUSCO_TSVTOGFF.out.gff.map { _meta, gff -> gff }.collect()
        )
        //ch_tree_data = ch_tree_data.mix(ORTHOSEQCOUNT.out.species_summary)
    }

    emit:
    orthofinder             = !val_skip_busco ? ch_orthofinder : channel.empty()        // channel: [ val(meta), [folder] ]
    tree_data               = !val_skip_busco ? ch_tree_data.flatten().collect().map { files -> files.toSorted { a, b -> a.name <=> b.name } } : channel.empty() // sort for deterministic behaviour
    quast_results           = QUAST.out.results                   // channel: [ val(meta), [tsv] ]
    busco_short_summaries   = !val_skip_busco ? BUSCO_GENOME.out.short_summaries_txt : channel.empty() // channel: [ val(meta), [txt] ]
    buscos_per_seqs         = !val_skip_busco ? IDEOGRAM_GENOME.out.busco_mappings.collect { _meta, table -> table}.map { tables -> tables.toSorted { a, b -> a.name <=> b.name } } : channel.empty() // channel: [ csv ]
    busco_batch_summaries   = !val_skip_busco ? GAWK.out.output : channel.empty()                     // channel: [ val(meta), path(tsv) ] — species-named batch summaries
    quast_tsv               = QUAST.out.tsv                       // channel: [ val(meta), path(tsv) ]

}
