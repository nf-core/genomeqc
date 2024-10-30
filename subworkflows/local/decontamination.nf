include { FCS_FCSADAPTOR               } from '../../modules/nf-core/fcsadaptor/main'
include { FCS_CLEANADAPTOR             } from '../../modules/local/fcs/cleanadaptor/main'
include { FCS_FCSGX                    } from '../../modules/nf-core/fcs/fcsgx/main'
include { FCS_CLEANGENOME              } from '../../modules/local/fcs/cleangenome/main'
include { TIARA_TIARA                  } from '../../modules/nf-core/tiara/tiara/main'



workflow DECONTAMINATION {

    take:
    // TODO nf-core: edit input (take) channels
    ch_fasta      // channel: [ val(meta), [ fasta ] ]

    main:

    ch_versions = Channel.empty()

    // Run module FCS_Adaptor find contamination

    FCS_FCSADAPTOR ( ch_fasta )
    ch_versions = ch_versions.mix(FCS_FCSADAPTOR.out.versions.first())

    //
    // Run Module Clean adaptor contamination
    //

    //ch_cleanadaptor_in = ch_fasta.join(FCS_FCSADAPTOR.out.adaptor_report, by: 0).view()

    ch_cleanadaptor_in = ch_fasta.join(FCS_FCSADAPTOR.out.adaptor_report, by: 0)
    .map { meta, fasta, adaptor_report -> tuple(meta, fasta, adaptor_report) }
    .view()

    FCS_CLEANADAPTOR (ch_cleanadaptor_in)
    //ch_versions = ch_versions.mix(FCS_CLEANADAPTOR.out.versions.first())


    FCS_FCSGX ( FCS_CLEANADAPTOR.out.cleaned,
    params.gxdb)
    ch_versions = ch_versions.mix(FCS_FCSGX.out.versions.first())



    ch_cleangenome_in = ch_fasta.join(FCS_FCSGX.out.fcs_gx_report, by: 0)
    .map { meta, fasta, fcs_gx_report -> tuple(meta, fasta, fcs_gx_report) }
    .view()

    FCS_CLEANGENOME (ch_cleangenome_in)
    //ch_versions = ch_versions.mix(FCSGX_CLEANGENOME.out.versions.first())

    TIARA_TIARA (FCS_CLEANGENOME.out.cleaned)
    //ch_versions = ch_versions.mix(TIARA_TIARA.out.versions.first())




    emit:
    // TODO nf-core: edit emitted channels
    //clean_fasta      = FCS_CLEANADAPTOR.out.clean_fasta           // channel: [ val(meta), [ clean_fasta ] ]
    //fcs_gx_report    = FCS_FCSGX.out.fcs_gx_report          // channel: [ val(meta), [ fcs_gx_report ] ]
    //adaptor_report   = FCS_FCSADAPTOR.out.adaptor_report          // channel: [ val(meta), [ adaptor_report ] ]
    //cleaned          = FCSGX_CLEANGENOME.out.cleaned            // channel: [val(meta),    [cleaned]

    versions = ch_versions                     // channel: [ versions.yml ]
}

