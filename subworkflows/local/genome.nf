
include { QUAST                               } from '../../modules/nf-core/quast/main'
include { BUSCO_BUSCO                         } from '../../modules/nf-core/busco/busco/main'

workflow GENOME {

    take:
    ch_fasta // channel: [ val(meta), [ fasta ] ]

    main:

    ch_versions = Channel.empty()

    QUAST ( 
        ch_fasta,
        [[],[]],
        [[],[]]
    )
    ch_versions = ch_versions.mix(QUAST.out.versions.first())

    BUSCO_BUSCO (
        ch_fasta,
        "genome", // hardcoded, other options ('proteins', 'transcriptome') make no sense
        params.busco_lineage,
        params.busco_lineages_path ?: [],
        params.busco_config ?: []
    )
    ch_versions = ch_versions.mix(BUSCO_BUSCO.out.versions.first())

    emit:
    quast = QUAST.out.tsv
    busco = BUSCO_BUSCO.out.batch_summary

    versions = ch_versions                     // channel: [ versions.yml ]
}

