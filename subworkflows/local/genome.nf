
include { QUAST                               } from '../../modules/nf-core/quast/main'
include { BUSCO_BUSCO                         } from '../../modules/nf-core/busco/busco/main'
include { GENOME_BUSCO_IDEOGRAM               } from '../../modules/local/genome_ideogram'

workflow GENOME {

    take:
    ch_fasta // channel: [ val(meta), [ fasta ] ]

    main:

    ch_versions   = Channel.empty()

    QUAST ( 
        ch_fasta,
        [[],[]],
        [[],[]]
    )
    ch_versions   = ch_versions.mix(QUAST.out.versions.first())

    BUSCO_BUSCO (
        ch_fasta,
        "genome", // hardcoded, other options ('proteins', 'transcriptome') make no sense
        params.busco_lineage,
        params.busco_lineages_path ?: [],
        params.busco_config ?: []
    )
    ch_versions   = ch_versions.mix(BUSCO_BUSCO.out.versions.first())
    ch_full_table = BUSCO_BUSCO.out.full_table

    // Combined ch_fasta and BUSCO output channel into a single channel for ideogram
    ch_input_ideo = ch_fasta
                  | combine(ch_full_table, by:0)


    GENOME_BUSCO_IDEOGRAM (
        ch_input_ideo
    )

    emit:
    quast_results         = QUAST.out.results                   // channel: [ val(meta), [tsv] ]
    busco_short_summaries = BUSCO_BUSCO.out.short_summaries_txt // channel: [ val(meta), [txt] ]

    versions = ch_versions                                      // channel: [ versions.yml ]
}

