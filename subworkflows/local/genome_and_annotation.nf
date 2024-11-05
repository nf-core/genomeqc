
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

    // TODO nf-core: substitute modules here for the modules of your subworkflow

    // Check GFF integrity
    ch_agat_gff = AGAT_CONVERTSPGXF2GXF(ch_gff).output_gff

    //
    // Run Quast
    //

    QUAST (
        ch_fasta,
        [[],[]],
        ch_agat_gff
    )
    ch_versions = ch_versions.mix(QUAST.out.versions.first())

    // For tree

    ch_tree_data = ch_tree_data.mix(QUAST.out.tsv.map { tuple -> tuple[1] })

    //
    // Run AGAT Spstatistics
    //

    AGAT_SPSTATISTICS (
        ch_agat_gff
    )
    ch_versions = ch_versions.mix(AGAT_SPSTATISTICS.out.versions.first())

    //
    // Run AGAT longest isoform
    //

//    LONGEST (
//        ch_ch_agat_gff
//    )
//    ch_versions = ch_versions.mix(LONGEST.out.versions.first())
//
//    //
//    // Run GFFREAD
//    //
//
//    ch_long_gff = LONGEST.out.longest_proteins
//    
    inputChannel = ch_agat_gff.combine(ch_fasta, by: 0)

    // Split the input channel into two channels
    gffChannel = inputChannel.map { tuple ->
        // Extracting the GFF path and ID
        [tuple[0], tuple[1]]
    }
    fnaChannel = inputChannel.map { tuple ->
        // Extracting only the FNA path
        [tuple[0], tuple[2]]
    }

    GFFREAD ( 
        fnaChannel,
        gffChannel
    )
    ch_versions = ch_versions.mix(GFFREAD.out.versions.first())

    //
    // MODULE: Run Orthofinder
    //

    ortho_ch = GFFREAD.out.longest.collect().map { it -> [[id:"orthofinder"], it] }

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
        "proteins", // Hard coded, it's the only possible option
        params.busco_lineage,
        params.busco_lineages_path ?: [],
        params.busco_config ?: []
    )
    ch_versions = ch_versions.mix(BUSCO_BUSCO.out.versions.first())

    // Prepare BUSCO output
    ch_busco_full_table = BUSCO_BUSCO.out.full_table
        .map { meta, full_tables -> 
            def genus = meta.id.split('_')[0]
            def speci = meta.id.split('_')[1]
            def genusspeci = "${genus}_${speci}"
            def lineages = full_tables.collect { it.toString().split('/')[-2].replaceAll('run_', '').replaceAll('_odb\\d+', '') }
            [genusspeci, lineages, full_tables]
        }
        //.view()

    // Prepare AGAT output
    ch_agat_gff = ch_agat_gff
        .map { meta, gff -> 
            def genus = meta.id.split('_')[0]
            def speci = meta.id.split('_')[1]
            def genusspeci = "${genus}_${speci}"
            [genusspeci, gff]
        }
        //.view()

    // Combine BUSCO and AGAT outputs
    ch_plot_input = ch_busco_full_table
        .join(ch_agat_gff)
        .flatMap { genusspeci, lineages, full_tables, gff ->
            lineages.withIndex().collect { lineage, index ->
                [genusspeci, lineage, full_tables[index], gff]
            }
        }

    //ch_plot_input.view()

    PLOT_BUSCO_IDEOGRAM ( ch_plot_input )//removed this temporarily:, ch_karyotype

    ch_tree_data = ch_tree_data.mix(BUSCO_BUSCO.out.batch_summary.collect { meta, file -> file })

    emit:
    orthofinder = ORTHOFINDER.out.orthofinder // channel: [ val(meta), [folder] ]
    //busco = BUSCO_BUSCO.out.batch_summary.collect { meta, file -> file }
    tree_data = ch_tree_data.flatten().collect()

    versions = ch_versions                     // channel: [ versions.yml ]
}

