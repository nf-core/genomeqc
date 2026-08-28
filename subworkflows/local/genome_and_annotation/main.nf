
include { AGAT_CONVERTSPGXF2GXF as AGAT_VALIDATE } from '../../../modules/nf-core/agat/convertspgxf2gxf'
include { AGAT_SPKEEPLONGESTISOFORM              } from '../../../modules/nf-core/agat/spkeeplongestisoform'
include { BUSCO_BUSCO as BUSCO_GENOME            } from '../../../modules/nf-core/busco/busco/main'
include { BUSCO_BUSCO as BUSCO_PROTEIN          } from '../../../modules/nf-core/busco/busco/main'
include { QUAST                                  } from '../../../modules/nf-core/quast/main'
include { AGAT_SPSTATISTICS                      } from '../../../modules/nf-core/agat/spstatistics/main'
include { IDEOGRAM_GENOMEANNOTATION              } from '../../../modules/local/ideogram/genomeannotation/main'
include { GFFREAD                                } from '../../../modules/nf-core/gffread/main'
include { GFFREAD as GFFREAD_VALIDATE            } from '../../../modules/nf-core/gffread/main'
include { ORTHOFINDER as ORTHOFINDER_V3          } from '../../../modules/nf-core/orthofinder/main'
include { ORTHOFINDERV2                          } from '../../../modules/local/orthofinderv2/main'
include { GENEOVERLAPS                           } from '../../../modules/local/geneoverlaps/main'
include { ORTHOSEQCOUNT                          } from '../../../modules/local/orthoseqcount'
include { GAWK as GAWK_GENO                      } from '../../../modules/nf-core/gawk/main'
include { GAWK as GAWK_PROT                      } from '../../../modules/nf-core/gawk/main'

workflow GENOME_AND_ANNOTATION {

    take:
    ch_fasta              // channel: [ val(meta), [ fasta ] ]
    ch_gxf                // channel: [ val(meta), [ gxf ] ]
    ch_busco_db           // channel: [ val(lineage) ]
    val_validation_tool   // val: GXF validation/standardisation tool - 'agat' or 'gffread'
    val_ortho_version     // val: OrthoFinder version - 'v2' or 'v3'
    val_skip_busco        // val: boolean - skip BUSCO (and everything downstream of it: ideogram)
    val_busco_lineage     // val: BUSCO lineage name (e.g. 'hymenoptera_odb10') or 'auto'
    val_busco_config      // val: path to a BUSCO config file, or []
    val_busco_clean       // val: boolean - clean up intermediate BUSCO files, or []

    main:
    ch_fasta.view { meta, _fasta -> "Running ${meta.id} on genome and annotation mode" }

    ch_versions  = channel.empty()

    // For tree plot
    ch_tree_data = channel.empty()

    //
    // MODULE: Run AGAT convertspgxf2gxf or GFFREAD validate
    //

    // Fix and standarize GXF
    if ( val_validation_tool == "agat" ) {
        AGAT_VALIDATE (
            ch_gxf
        )
        ch_gxf_agat  = AGAT_VALIDATE.out.output_gff
    } else if ( val_validation_tool == "gffread" ) {
        GFFREAD_VALIDATE (
            ch_gxf,
            []
        )
        ch_gxf_agat  = GFFREAD_VALIDATE.out.gffread_gff
    }

    //
    // MODULE: Run AGAT longest isoform
    //


    AGAT_SPKEEPLONGESTISOFORM (
        ch_gxf_agat,
        []

    )

    // Get longest isoform from gff
    ch_gxf_long  = AGAT_SPKEEPLONGESTISOFORM.out.gff


    //
    // Prepare input multichannel
    //

    // Combine inputs (fasta, gff from AGAT (unfiltered) and gff from AGAT_SPKEEPLONGESTISOFORM (filtered)))
    // into a single multichannel so that they are in sync
    ch_input     = ch_fasta
        | combine(ch_gxf_agat, by:0) // by:0 | Only combine when both channels share the same id
        | combine(ch_gxf_long, by:0)
        | multiMap {
            meta, fasta, gxf_unfilt, gxf_filt -> // "null" probably not necessary
                fasta      : fasta      ? tuple( meta, file(fasta)      ) : null // channel: [ val(meta), [ fasta ] ]
                gxf_unfilt : gxf_unfilt ? tuple( meta, file(gxf_unfilt) ) : null // channel: [ val(meta), [ gxf ] ], unfiltered
                gxf_filt   : gxf_filt   ? tuple( meta, file(gxf_filt)   ) : null // channel: [ val(meta), [ gxf ] ], filtered for longest isoform
        }

    //
    // Run AGAT Spstatistics
    //

    AGAT_SPSTATISTICS (
        ch_input.gxf_unfilt
    )

    //
    // MODULE: Run gene overlap module
    //

    GENEOVERLAPS {
        ch_input.gxf_filt
    }
    ch_tree_data = ch_tree_data.mix(GENEOVERLAPS.out.overlap_counts.collect { _meta, file -> file })

    //
    // MODULE: Run Quast
    //

    QUAST (
        ch_input.fasta,
        [[],[]],
        ch_input.gxf_unfilt
    )

    // For tree

    ch_tree_data = ch_tree_data.mix(QUAST.out.tsv.map { tuple -> tuple[1] })

    //
    // MODULE: Run GFFREAD
    //

    GFFREAD (
        ch_input.gxf_filt,
        ch_input.fasta.map { _meta, fasta -> fasta}
    )

    //
    // MODULE: Run Orthofinder
    //

    // Prepare orthofinder input channel
    ortho_ch     = GFFREAD.out.gffread_fasta
        | map { _meta, fasta ->
            fasta // We only need the fastas
        }
        | collect // Collect all fasta in a single tuple
        | filter { fastas ->
            fastas.size() >= 3 // Ensure we have at least 3 genomes for orthofinder, otherwise it won't run
        }
        | map { fastas ->
            [[id:'orthofinder', mode:'genome_anno'], fastas.toSorted { a, b -> a.name <=> b.name }] // sort for deterministic behaviour
        }

    // Run orthofinder
    if (val_ortho_version == 'v3') {
        ORTHOFINDER_V3 (
            ortho_ch,
            [[],[]]
        )
        ch_orthofinder = ORTHOFINDER_V3.out.orthofinder
    } else if (val_ortho_version == 'v2' ) {
        ORTHOFINDERV2 (
            ortho_ch
        )
        ch_orthofinder = ORTHOFINDERV2.out.orthofinder
    }

    //
    // MODULE: Run ORTHOSEQCOUNT
    //

    ORTHOSEQCOUNT (
        ch_orthofinder.map { _meta, folder ->
            file("${folder}/Orthogroups/Orthogroups.tsv")
        },
        AGAT_SPKEEPLONGESTISOFORM.out.gff.map { _meta, gff -> gff }.collect()
    )
    //ch_tree_data = ch_tree_data.mix(ORTHOSEQCOUNT.out.species_summary)

    //
    // MODULE: Run BUSCO for genome annotation
    //

    if(!val_skip_busco) {
        BUSCO_GENOME (
            ch_fasta,
            'genome',
            val_busco_lineage,
            ch_busco_db,
            val_busco_config ?: [],
            val_busco_clean ?: []
        )

        //
        // MODULE: Run BUSCO for proteins
        //

        BUSCO_PROTEIN (
            GFFREAD.out.gffread_fasta,
            'proteins',
            val_busco_lineage,
            ch_busco_db,
            val_busco_config ?: [],
            val_busco_clean ?: []
        )

        //
        // GAWK
        //
        // Use GAWK to change ID from file name to meta.id
        // For BUSCO genome
        GAWK_GENO (
            BUSCO_GENOME.out.batch_summary,
            [],
            false
        )

        // For BUSCO protein
        GAWK_PROT (
            BUSCO_PROTEIN.out.batch_summary,
            [],
            false
        )

        //
        // Plot BUSCO ideogram
        //

        // Prepare BUSCO output
        ch_busco_full_table = BUSCO_PROTEIN.out.full_table
                            | map { meta, full_tables ->
                                def lineages = full_tables.toString().split('/')[-2].replaceAll('run_', '').replaceAll('_odb\\d+', '')
                                [meta.id, lineages, full_tables]
                            }
                            | groupTuple(by: 0)
                            | map { id, lineages, full_tables ->
                                [id, lineages.flatten(), full_tables.flatten()]
                            }

        // Add genome to channel
        fnaChannel_busco    = ch_input.fasta
                            | map { meta, fasta ->
                                [meta.id, fasta]
                            }

        // Prepare GXF channel of ideogram
        ch_gxf_busco        = ch_input.gxf_filt
                            | map { meta, gxf ->
                                [meta.id, gxf]
                            }

        // Combine BUSCO, AGAT, and genome outputs
        ch_plot_input       = ch_busco_full_table
                            | join(fnaChannel_busco)
                            | join(ch_gxf_busco)
                            | flatMap { genusspeci, lineages, full_tables, fasta, gxf ->
                                lineages.withIndex().collect { lineage, index ->
                                    [ [id: genusspeci, lineage: lineage], full_tables[index], fasta, gxf ]
                                }
                            }

        IDEOGRAM_GENOMEANNOTATION ( ch_plot_input )
    }

    emit:
    orthofinder                = ch_orthofinder         // channel: [ val(meta), [folder] ]
    tree_data                  = ch_tree_data.flatten().collect().map { files -> files.toSorted { a, b -> a.name <=> b.name } } // sort for deterministic behaviour
    quast_results              = QUAST.out.results                   // channel: [ val(meta), [tsv] ]
    busco_short_summaries_geno = !val_skip_busco ? GAWK_GENO.out.output : channel.empty()
    busco_short_summaries_prot = !val_skip_busco ? GAWK_PROT.out.output : channel.empty()
    quast_tsv                  = QUAST.out.tsv                       // channel: [ val(meta), path(tsv) ]
    agat_stats                 = AGAT_SPSTATISTICS.out.stats_txt     // channel: [ val(meta), path(txt) ]
    ortho_seq_count            = ORTHOSEQCOUNT.out.species_summary // channel: [ path(tsv) ]
    buscos_per_seqs            = !val_skip_busco ? IDEOGRAM_GENOMEANNOTATION.out.busco_mappings.collect { _meta, table -> table}.map { tables -> tables.toSorted { a, b -> a.name <=> b.name } } : channel.empty() // channel: [ val(meta), [csv] ]
    versions                   = ch_versions                   // channel: [ versions.yml ]
}
