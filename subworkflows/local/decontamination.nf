// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { FCS_FCSADAPTOR     } from '../modules/nf-core/fcs/fcsadaptor/main'
include { FCS_CLEANADAPTOR   } from '../modules/local/fcs/cleanadaptor/main'
include { FCS_FCSGX          } from '../modules/nf-core/fcs/fcsgx/main'
include { FCSGX_CLEANGENOME  } from '../modules/local/fcs/cleangenome/main'
include { TIARA_TIARA        } from '../modules/nf-core/tiara/tiara/main'



workflow DECONTAMINATION {

    take:
    // TODO nf-core: edit input (take) channels
    ch_assembly      // channel: [ val(meta), [ fasta ] ]

    main:

    ch_versions = Channel.empty()

    // TODO nf-core: substitute modules here for the modules of your subworkflow

    FCS_FCSADAPTOR ( ch_assembly )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    FCS_FCSGX ( SAMTOOLS_SORT.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    emit:
    // TODO nf-core: edit emitted channels
    bam      = SAMTOOLS_SORT.out.bam           // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    csi      = SAMTOOLS_INDEX.out.csi          // channel: [ val(meta), [ csi ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

