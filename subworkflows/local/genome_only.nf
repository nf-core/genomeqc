include { QUAST                               } from '../../modules/nf-core/quast/main'
include { BUSCO_BUSCO                         } from '../../modules/nf-core/busco/busco/main'
include { GENOME_ONLY_BUSCO_IDEOGRAM          } from '../../modules/local/genome_only_busco_ideogram'
include { ORTHOFINDER                         } from '../../modules/nf-core/orthofinder/main'
include { BUSCO_TSV_TO_GFF                    } from '../../modules/local/busco_tsv_to_gff/main'
include { ORTHOLOGOUS_CHROMOSOMES             } from '../../modules/local/orthologous_chromosomes'
include { GAWK                                } from '../../modules/nf-core/gawk/main'

workflow GENOME_ONLY {

    take:
    ch_fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_fasta.view { "Running ${it[0]} on genome only mode"}

    ch_versions   = Channel.empty()

    // For tree plot
    ch_tree_data = Channel.empty()

    //
    // MODULE: Run Quast
    //

    QUAST (
        ch_fasta,
        [[],[]],
        [[],[]]
    )
    ch_versions   = ch_versions.mix(QUAST.out.versions.first())
    ch_tree_data = ch_tree_data.mix(QUAST.out.tsv.map { tuple -> tuple[1] })

    BUSCO_BUSCO (
        ch_fasta,
        "genome", // hardcoded, other options ('proteins', 'transcriptome') make no sense
        params.busco_lineage,
        params.busco_lineages_path ?: [],
        params.busco_config ?: [],
        params.busco_clean ?: []
    )
    ch_versions   = ch_versions.mix(BUSCO_BUSCO.out.versions.first())
    //ch_tree_data  = ch_tree_data.mix(BUSCO_BUSCO.out.batch_summary.collect { meta, file -> file })

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
    // BUSCO Ideogram
    //

    ch_full_table = BUSCO_BUSCO.out.full_table

    // Combined ch_fasta and BUSCO output channel into a single channel for ideogram
    ch_input_ideo = ch_fasta
                  | combine(ch_full_table, by:0)


    GENOME_ONLY_BUSCO_IDEOGRAM (
        ch_input_ideo
    )
    ch_versions   = ch_versions.mix(GENOME_ONLY_BUSCO_IDEOGRAM.out.versions.first())

    //
    // Orthofinder
    //

    // Prepare data
    ch_busco_proteins = BUSCO_BUSCO.out.single_copy_faa
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
                                [[ id: 'orthofinder_run', mode: 'genome' ], file_paths]
                      }

    //Run orthofinder
    ORTHOFINDER (
        ch_busco_proteins,
        [[],[]]
    )
    ch_versions  = ch_versions.mix(ORTHOFINDER.out.versions)

    // Transform tsv to gff for orthologous chromosomes module
    BUSCO_TSV_TO_GFF (
        BUSCO_BUSCO.out.busco_dir
    )
    //ch_versions  = ch_versions.mix(BUSCO_TSV_TO_GFF.out.versions)

    //
    // MODULE: Run ORTHOLOGOUS_CHROMOSOMES
    //

    ORTHOLOGOUS_CHROMOSOMES (
        ORTHOFINDER.out.orthofinder.map { meta, folder ->
            file("${folder}/Orthogroups/Orthogroups.tsv")
        },
        BUSCO_TSV_TO_GFF.out.gff.map { meta, gff -> gff }.collect()
    )
    ch_versions  = ch_versions.mix(ORTHOLOGOUS_CHROMOSOMES.out.versions)
    ch_tree_data = ch_tree_data.mix(ORTHOLOGOUS_CHROMOSOMES.out.species_summary)

    emit:
    orthofinder             = ORTHOFINDER.out.orthofinder         // channel: [ val(meta), [folder] ]
    tree_data               = ch_tree_data.flatten().collect()
    quast_results           = QUAST.out.results                   // channel: [ val(meta), [tsv] ]
    busco_short_summaries   = BUSCO_BUSCO.out.short_summaries_txt // channel: [ val(meta), [txt] ]
    buscos_per_seqs         = GENOME_ONLY_BUSCO_IDEOGRAM.out.busco_mappings.collect { meta, table -> table} // channel: [ csv ]

    versions = ch_versions                                      // channel: [ versions.yml ]
}
