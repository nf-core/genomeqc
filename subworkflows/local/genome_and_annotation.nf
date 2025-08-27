
include { AGAT_CONVERTSPGXF2GXF               } from '../../modules/nf-core/agat/convertspgxf2gxf'
include { AGAT_SPKEEPLONGESTISOFORM           } from '../../modules/nf-core/agat/spkeeplongestisoform'
include { BUSCO_BUSCO                         } from '../../modules/nf-core/busco/busco/main'
include { QUAST                               } from '../../modules/nf-core/quast/main'
include { AGAT_SPSTATISTICS                   } from '../../modules/nf-core/agat/spstatistics/main'
include { GENOME_ANNOTATION_BUSCO_IDEOGRAM    } from '../../modules/local/genome_annotation_busco_ideogram'
include { GFFREAD                             } from '../../modules/nf-core/gffread/main'
include { GFFREAD as GFFREAD_VALIDATE         } from '../../modules/nf-core/gffread/main'
include { ORTHOFINDER                         } from '../../modules/nf-core/orthofinder/main'
include { FASTAVALIDATOR                      } from '../../modules/nf-core/fastavalidator/main'
include { GENE_OVERLAPS                       } from '../../modules/local/gene_overlaps'
include { ORTHOLOGOUS_CHROMOSOMES             } from '../../modules/local/orthologous_chromosomes'
include { GAWK                                } from '../../modules/nf-core/gawk/main'

workflow GENOME_AND_ANNOTATION {

    take:
    ch_fasta // channel: [ val(meta), [ fasta ] ]
    ch_gxf   // channel: [ val(meta), [ gxf ] ]

    main:
    ch_fasta.view { "Running ${it[0]} on genome and annotation mode"}

    ch_versions  = Channel.empty()

    // For tree plot
    ch_tree_data = Channel.empty()

    //
    // MODULE: run AGAT convertspgxf2gxf
    //

    // Fix and standarize GXF
    //AGAT_CONVERTSPGXF2GXF(
    //    ch_gxf
    //)
    //ch_gxf_agat  = AGAT_CONVERTSPGXF2GXF.out.output_gff
    //ch_versions  = ch_versions.mix(AGAT_CONVERTSPGXF2GXF.out.versions.first())

    GFFREAD_VALIDATE (
        ch_gxf,
        []
    )
    ch_gxf_agat  = GFFREAD_VALIDATE.out.gffread_gff
    

    //
    // MODULE: Run AGAT longest isoform
    //


    AGAT_SPKEEPLONGESTISOFORM (
        ch_gxf_agat,
        []

    )
    ch_versions  = ch_versions.mix(AGAT_SPKEEPLONGESTISOFORM.out.versions.first())

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
    ch_versions  = ch_versions.mix(AGAT_SPSTATISTICS.out.versions.first())

    //
    // MODULE: Run gene overlap module
    //

    GENE_OVERLAPS {
        ch_input.gxf_filt
    }
    ch_versions  = ch_versions.mix(GENE_OVERLAPS.out.versions.first())
    ch_tree_data = ch_tree_data.mix(GENE_OVERLAPS.out.overlap_counts.collect { meta, file -> file })

    //
    // MODULE: Run Quast
    //

    QUAST (
        ch_input.fasta,
        [[],[]],
        ch_input.gxf_unfilt
    )
    ch_versions  = ch_versions.mix(QUAST.out.versions.first())

    // For tree

    ch_tree_data = ch_tree_data.mix(QUAST.out.tsv.map { tuple -> tuple[1] })

    //
    // MODULE: Run GFFREAD
    //

    GFFREAD (
        ch_input.gxf_filt,
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
    ch_versions  = ch_versions.mix(FASTAVALIDATOR.out.versions.first())

    //
    // MODULE: Run Orthofinder
    //

    // Prepare orthofinder input channel
    ortho_ch     = GFFREAD.out.gffread_fasta
        | map { meta, fasta ->
            fasta // We only need the fastas
        }
        | collect // Collect all fasta in a single tuple
        | filter { fastas ->
            fastas.size() >= 3 // Ensure we have at least 4 genomes for orthofinder, otherwise it won't run
        }
        | map { fastas ->
            [[id:'orthofinder', mode:'genome_anno'], fastas]
        }

    // Run orthofinder
    ORTHOFINDER (
        ortho_ch,
        [[],[]]
    )
    ch_versions  = ch_versions.mix(ORTHOFINDER.out.versions)
    ORTHOFINDER.out.orthofinder.view()

    //
    // MODULE: Run ORTHOLOGOUS_CHROMOSOMES
    //

    ORTHOLOGOUS_CHROMOSOMES (
        ORTHOFINDER.out.orthofinder.map { meta, folder ->
            file("${folder}/Orthogroups/Orthogroups.tsv")
        },
        AGAT_SPKEEPLONGESTISOFORM.out.gff.map { meta, gff -> gff }.collect()
    )
    ch_versions  = ch_versions.mix(ORTHOLOGOUS_CHROMOSOMES.out.versions)
    ch_tree_data = ch_tree_data.mix(ORTHOLOGOUS_CHROMOSOMES.out.species_summary)

    //
    // MODULE: Run BUSCO
    //

    BUSCO_BUSCO (
        GFFREAD.out.gffread_fasta,
        params.busco_mode,
        params.busco_lineage,
        params.busco_lineages_path ?: [],
        params.busco_config ?: [],
        params.busco_clean ?: []
    )
    ch_versions  = ch_versions.mix(BUSCO_BUSCO.out.versions.first())

    //
    // GAWK
    //
    // Use GAWK to change ID from file name to meta.id

    GAWK (
        BUSCO_BUSCO.out.batch_summary,
        [],
        false
    )
    ch_tree_data  = ch_tree_data.mix(GAWK.out.output.collect { meta, file -> file })

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
                                [genusspeci, lineage, full_tables[index], fasta, gxf]
                            }
                        }

    GENOME_ANNOTATION_BUSCO_IDEOGRAM ( ch_plot_input )
    ch_versions         = ch_versions.mix(GENOME_ANNOTATION_BUSCO_IDEOGRAM.out.versions.first())

    //ch_tree_data        = ch_tree_data.mix(BUSCO_BUSCO.out.batch_summary.collect { meta, file -> file })

    emit:
    orthofinder           = ORTHOFINDER.out.orthofinder         // channel: [ val(meta), [folder] ]
    tree_data             = ch_tree_data.flatten().collect()
    quast_results         = QUAST.out.results                   // channel: [ val(meta), [tsv] ]
    busco_short_summaries = BUSCO_BUSCO.out.short_summaries_txt // channel: [ val(meta), [txt] ]
    orthologous_chromosomes = ORTHOLOGOUS_CHROMOSOMES.out.species_summary // channel: [ path(tsv) ]
    buscos_per_seqs       = GENOME_ANNOTATION_BUSCO_IDEOGRAM.out.busco_mappings.collect { meta, table -> table} // channel: [ val(meta), [csv] ]

    versions              = ch_versions                   // channel: [ versions.yml ]
}
