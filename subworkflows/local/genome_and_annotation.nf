
include { AGAT_CONVERTSPGXF2GXF               } from '../../modules/nf-core/agat/convertspgxf2gxf'
include { LONGEST                             } from '../../modules/local/longest'
include { BUSCO_BUSCO                         } from '../../modules/nf-core/busco/busco/main'
include { QUAST                               } from '../../modules/nf-core/quast/main'
include { AGAT_SPSTATISTICS                   } from '../../modules/nf-core/agat/spstatistics/main'
include { GENOME_ANNOTATION_BUSCO_IDEOGRAM    } from '../../modules/local/genome_annotation_busco_ideogram'
include { GFFREAD                             } from '../../modules/nf-core/gffread/main'
include { ORTHOFINDER                         } from '../../modules/nf-core/orthofinder/main'
include { FASTAVALIDATOR                      } from '../../modules/nf-core/fastavalidator/main'
include { GENE_OVERLAPS                       } from '../../modules/local/gene_overlaps'

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
                     meta, fasta, gff_unfilt, gff_filt -> // "null" probably not necessary
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
    // MODULE: Run gene overlap module
    //

    GENE_OVERLAPS {
        ch_input.gff_filt
    }
    ch_tree_data = ch_tree_data.mix(GENE_OVERLAPS.out.overlap_counts.collect { meta, file -> file })

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
    // MODULE: Run GFFREAD
    //

    GFFREAD (
        ch_input.gff_filt,
        ch_input.fasta.map { meta, fasta -> fasta}
    )
    ch_versions  = ch_versions.mix(GFFREAD.out.versions.first())

    //
    // MODULE: Run fasta validator
    //

    // Shoud we keep this?
    FASTAVALIDATOR(
        GFFREAD.out.gffread_fasta
    )

    //
    // MODULE: Run Orthofinder
    //

    // Prepare orthofinder input channel
    ortho_ch     = GFFREAD.out.gffread_fasta
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
    // MODULE: Run BUSCO
    //

    BUSCO_BUSCO (
        GFFREAD.out.gffread_fasta,
        "proteins", // hardcoded
        params.busco_lineage,
        params.busco_lineages_path ?: [],
        params.busco_config ?: []
    )
    ch_versions  = ch_versions.mix(BUSCO_BUSCO.out.versions.first())

    //
    // Plot BUSCO ideogram
    //

    // Prepare BUSCO output
    ch_busco_full_table = BUSCO_BUSCO.out.full_table
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

    GENOME_ANNOTATION_BUSCO_IDEOGRAM ( ch_plot_input )

    ch_tree_data        = ch_tree_data.mix(BUSCO_BUSCO.out.batch_summary.collect { meta, file -> file })

    emit:
    orthofinder           = ORTHOFINDER.out.orthofinder         // channel: [ val(meta), [folder] ]
    tree_data             = ch_tree_data.flatten().collect()
    quast_results         = QUAST.out.results                   // channel: [ val(meta), [tsv] ]
    busco_short_summaries = BUSCO_BUSCO.out.short_summaries_txt // channel: [ val(meta), [txt] ]

    versions              = ch_versions                         // channel: [ versions.yml ]
}
