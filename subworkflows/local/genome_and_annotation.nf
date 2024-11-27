
include { AGAT_CONVERTSPGXF2GXF               } from '../../modules/nf-core/agat/convertspgxf2gxf'
include { LONGEST                             } from '../../modules/local/longest'
include { QUAST                               } from '../../modules/nf-core/quast/main'
include { AGAT_SPSTATISTICS                   } from '../../modules/nf-core/agat/spstatistics/main'
include { PLOT_BUSCO_IDEOGRAM                 } from '../../modules/local/plot_busco_ideogram'
include { ORTHOFINDER                         } from '../../modules/nf-core/orthofinder/main'
include { FASTAVALIDATOR                      } from '../../modules/nf-core/fastavalidator/main'

include { FASTA_GXF_BUSCO_PLOT                } from '../../subworkflows/nf-core/fasta_gxf_busco_plot/main'

workflow GENOME_AND_ANNOTATION {

    take:
    ch_fasta // channel: [ val(meta), [ fasta ] ]
    ch_gff   // channel: [ val(meta), [ gff ] ]

    main:

    ch_versions  = Channel.empty()

    // For tree plot
    ch_tree_data = Channel.empty()

    //
    // MODULE: run AGAT convertspgxf2gxf
    //

    // Fix and standarize GFF
    ch_gff_agat  = AGAT_CONVERTSPGXF2GXF(ch_gff).output_gff

    //
    // MODULE: Run AGAT longest isoform
    //

    LONGEST (
        ch_gff_agat
    )
    ch_versions  = ch_versions.mix(LONGEST.out.versions.first())

    // Get longest isoform from gff
    ch_gff_long  = LONGEST.out.longest_proteins

    //
    // Prepare input multichannel
    //

    // Combine inputs (fasta, gff from AGAT (unfiltered) and gff from LONGEST (fitlered)))
    // into a single multichannel so that they are in sync
    ch_input     = ch_fasta
                 | combine(ch_gff_agat, by:0) // by:0 | Only combine when both channels share the same id
                 | combine(ch_gff_long, by:0)
                 | multiMap {
                     meta, fasta, gff_unfilt, gff_filt ->
                         fasta      : fasta      ? tuple( meta, file(fasta)      ) : null // channel: [ val(meta), [ fasta ] ]
                         gff_unfilt : gff_unfilt ? tuple( meta, file(gff_unfilt) ) : null // channel: [ val(meta), [ gff ] ], unfiltered
                         gff_filt   : gff_filt   ? tuple( meta, file(gff_filt)   ) : null // channel: [ val(meta), [ gff ] ], filtered for longest isoform
                 }

    //
    // Run AGAT Spstatistics
    //

    AGAT_SPSTATISTICS (
        ch_input.gff_unfilt
    )
    ch_versions  = ch_versions.mix(AGAT_SPSTATISTICS.out.versions.first())

    //
    // MODULE: Run Quast
    //

    QUAST (
        ch_input.fasta,
        [[],[]],
        ch_input.gff_unfilt
    )
    ch_versions  = ch_versions.mix(QUAST.out.versions.first())

    // For tree

    ch_tree_data = ch_tree_data.mix(QUAST.out.tsv.map { tuple -> tuple[1] })

    //
    // SUBWORKFLOW: FASTA_GXF_BUSCO_PLOT
    //
    FASTA_GXF_BUSCO_PLOT (
        ch_input.fasta,
        ch_input.gff_filt,
        'genome', // mode
        [ params.busco_lineage ], // List expected
        params.busco_lineages_path,
        params.busco_config
    )

    ch_versions  = ch_versions.mix(FASTA_GXF_BUSCO_PLOT.out.versions)

    //
    // MODULE: Run fasta validator
    //

    // Shoud we keep this?
    FASTAVALIDATOR(
        FASTA_GXF_BUSCO_PLOT.out.proteins
    )

    //
    // MODULE: Run Orthofinder
    //

    // Prepare orthofinder input channel
    ortho_ch     = FASTA_GXF_BUSCO_PLOT.out.proteins
                 | map { meta, fasta ->
                     fasta // We only need the fastas
                 }
                 | collect // Collect all fasta in a single tuple
                 | map { fastas ->
                     [[id:"orthofinder"], fastas] 
                 }

    // Run orthofinder
    ORTHOFINDER (
        ortho_ch,
        [[],[]]
    )
    ch_versions  = ch_versions.mix(ORTHOFINDER.out.versions)

    //
    // Plot BUSCO ideogram
    //

    // Prepare BUSCO output
    FASTA_GXF_BUSCO_PLOT.out.annotation_full_table.map { meta, full_tables -> full_tables }.view()
    
    //BUSCO_BUSCO.out.full_table.map { meta, full_tables -> full_tables.collect { it -> it.toString().split('/')[-2] } }.view()

    //ch_busco_full_table = BUSCO_BUSCO.out.full_table
    //                        | map { meta, full_tables ->
    //                            def lineages = full_tables.collect { it -> it.toString().split('/')[-2].replaceAll('run_', '').replaceAll('_odb\\d+', '') }
    //                            [meta.id, lineages, full_tables]
    //                        }

    ch_busco_full_table = FASTA_GXF_BUSCO_PLOT.out.annotation_full_table
                            .map { meta, full_tables ->
                                def lineages = full_tables.toString().split('/')[-2].replaceAll('run_', '').replaceAll('_odb\\d+', '')
                                [meta.id, lineages, full_tables]
                            }
                            .groupTuple(by: 0)
                            .map { id, lineages, full_tables ->
                                [id, lineages.flatten(), full_tables.flatten()]
                            }

    // Add genome to channel
    fnaChannel_busco    = ch_input.fasta
                        | map { meta, fasta ->
                            [meta.id, fasta]
                        }

    // Prepare GFF channel of ideogram
    ch_gff_busco        = ch_input.gff_filt
                        | map { meta, gff ->
                            [meta.id, gff]
                        }

    // Combine BUSCO, AGAT, and genome outputs
    ch_plot_input       = ch_busco_full_table
                        | join(fnaChannel_busco)
                        | join(ch_gff_busco)
                        | flatMap { genusspeci, lineages, full_tables, fasta, gff ->
                            lineages.withIndex().collect { lineage, index ->
                                [genusspeci, lineage, full_tables[index], fasta, gff]
                            }
                        }

    PLOT_BUSCO_IDEOGRAM ( ch_plot_input )//removed this temporarily:, ch_karyotype

    ch_tree_data        = ch_tree_data.mix(FASTA_GXF_BUSCO_PLOT.out.annotation_batch_summary.collect { meta, file -> file })

    emit:
    orthofinder         = ORTHOFINDER.out.orthofinder // channel: [ val(meta), [folder] ]
    //busco = BUSCO_BUSCO.out.batch_summary.collect { meta, file -> file }

    tree_data           = ch_tree_data.flatten().collect()
    busco_mq            = FASTA_GXF_BUSCO_PLOT.out.annotation_short_summaries_txt.map { meta, file -> file } 
    quast_mq            = QUAST.out.results.map { meta, file -> file }

    versions            = ch_versions // channel: [ versions.yml ]
}
