
include { AGAT_CONVERTSPGXF2GXF               } from '../../modules/nf-core/agat/convertspgxf2gxf'
include { LONGEST                             } from '../../modules/local/longest'
include { BUSCO_BUSCO                         } from '../../modules/nf-core/busco/busco/main'
include { QUAST                               } from '../../modules/nf-core/quast/main'
include { AGAT_SPSTATISTICS                   } from '../../modules/nf-core/agat/spstatistics/main'
include { PLOT_BUSCO_IDEOGRAM                 } from '../../modules/local/plot_busco_ideogram.nf'
//include { GFFREAD                             } from '../../modules/nf-core/gffread/main'
include { GFFREAD                             } from '../../modules/local/gffread'
include { ORTHOFINDER                         } from '../../modules/nf-core/orthofinder/main'

workflow GENOME_AND_ANNOTATION {

    take:
    ch_fasta // channel: [ val(meta), [ fasta ] ]
    ch_gff // channel: [ val(meta), [ gff ] ]

    main:

    ch_versions = Channel.empty()
 
    // For tree plot
    ch_tree_data = Channel.empty()

    // Fix and standarize GFF
    ch_gff_agat = AGAT_CONVERTSPGXF2GXF(ch_gff).output_gff

    // Combine inputs into a single multichannel
    ch_fasta
        | combine(ch_gff_agat, by:0) // by:0 | Only combine when both channels share the same id
        | multiMap {
            meta, fasta, gff ->
                fasta : fasta ? tuple( meta, file(fasta) ) : null
                gff   : gff   ? tuple( meta, file(gff) )   : null
        }
        | set { ch_input }


    //
    // Run Quast
    //

    QUAST (
        ch_input.fasta,
        [[],[]],
        ch_input.gff
    )
    ch_versions = ch_versions.mix(QUAST.out.versions.first())

    // For tree

    ch_tree_data = ch_tree_data.mix(QUAST.out.tsv.map { tuple -> tuple[1] })

    //
    // Run GFFREAD
    //

    GFFREAD ( 
        ch_input.fasta,
        ch_input.gff
    )
    ch_versions = ch_versions.mix(GFFREAD.out.versions.first())

    //
    // MODULE: Run Orthofinder
    //

    ortho_ch = GFFREAD.out.longest.collect().map { fastas -> [[id:"orthofinder"], fastas] }
    
    ORTHOFINDER (
        ortho_ch,
        [[],[]]
    )
    ch_versions = ch_versions.mix(ORTHOFINDER.out.versions)

    //
    // MODULE: Run BUSCO
    //

    BUSCO_BUSCO (  
        GFFREAD.out.proteins_busco, 
        "proteins", // hard coded
        params.busco_lineage,
        params.busco_lineages_path ?: [],
        params.busco_config ?: []
    )
    ch_versions = ch_versions.mix(BUSCO_BUSCO.out.versions.first())

    // Prepare BUSCO output
    ch_busco_full_table = BUSCO_BUSCO.out.full_table
        .map { meta, full_tables -> 
            def lineages = full_tables.collect { it.toString().split('/')[-2].replaceAll('run_', '').replaceAll('_odb\\d+', '') }
            [meta.id, lineages, full_tables]
        }

    // Add genome to channel
    fnaChannel_busco = fnaChannel
        .map { meta, fasta -> 
            [meta.id, fasta]
        }

    // Prepare AGAT output
    ch_agat_gff_busco = ch_agat_gff
        .map { meta, gff -> 
            [meta.id, gff]
        }

    // Combine BUSCO, AGAT, and genome outputs
    ch_plot_input = ch_busco_full_table
        .join(fnaChannel_busco)
        .join(ch_agat_gff_busco)
        .flatMap { genusspeci, lineages, full_tables, fasta, gff ->
            lineages.withIndex().collect { lineage, index ->
                [genusspeci, lineage, full_tables[index], fasta, gff]
            }
        }

    PLOT_BUSCO_IDEOGRAM ( ch_plot_input )//removed this temporarily:, ch_karyotype

    ch_tree_data = ch_tree_data.mix(BUSCO_BUSCO.out.batch_summary.collect { meta, file -> file })

    //
    // Run AGAT Spstatistics
    //

    AGAT_SPSTATISTICS (
        ch_input.gff
    )
    ch_versions = ch_versions.mix(AGAT_SPSTATISTICS.out.versions.first())


    emit:
    orthofinder = ORTHOFINDER.out.orthofinder // channel: [ val(meta), [folder] ]
    //busco = BUSCO_BUSCO.out.batch_summary.collect { meta, file -> file }
    tree_data = ch_tree_data.flatten().collect()

    versions = ch_versions                     // channel: [ versions.yml ]
}

