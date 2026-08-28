include { FCS_FCSADAPTOR                          } from '../../../modules/nf-core/fcs/fcsadaptor/main'
include { FCSGX_CLEANGENOME as FCSGX_CLEANADAPTOR } from '../../../modules/nf-core/fcsgx/cleangenome/main'
include { FCSGX_FETCHDB                           } from '../../../modules/nf-core/fcsgx/fetchdb/main'
include { FCSGX_RUNGX                             } from '../../../modules/nf-core/fcsgx/rungx/main'
include { FCSGX_CLEANGENOME                       } from '../../../modules/nf-core/fcsgx/cleangenome/main'
include { TIARA_TIARA as TIARA_RAW                } from '../../../modules/nf-core/tiara/tiara/main'
include { TIARA_TIARA as TIARA_CLEANED            } from '../../../modules/nf-core/tiara/tiara/main'



workflow DECONTAMINATION {

    take:
    ch_fasta      // channel: [ val(meta), [ fasta ] ]
    ch_ramdisk
    ch_gxdb_local // channel: val(gxdb)
    ch_gxdb_manifest // channel: val(gxdb)

    main:

    // Run module FCS_Adaptor find contamination

    FCS_FCSADAPTOR (
        ch_fasta
    )

    //
    // Run Module Clean adaptor contamination
    //

    ch_cleanadaptor_in = ch_fasta
                       | join(FCS_FCSADAPTOR.out.adaptor_report, by: 0)
                       | map { meta, fasta, adaptor_report -> tuple(meta, fasta, adaptor_report) }

    FCSGX_CLEANADAPTOR (ch_cleanadaptor_in)

    // Run Module fetch database
    FCSGX_FETCHDB(
        channel.value(ch_gxdb_manifest).filter { manifest -> manifest != [] } // If there's no manifest, use empty channel (won't run)
    )
    ch_gxdb            =  ch_gxdb_local != [] ? channel.value(ch_gxdb_local) : FCSGX_FETCHDB.out.database // If there's no local database, use the fetched one

    // Run Module fcs gx

    // Prepare input channel for fcs gx
    ch_fcsgx = FCSGX_CLEANADAPTOR.out.cleaned
             | map { meta, fasta -> tuple( meta, meta.taxid, fasta) }

    FCSGX_RUNGX (
        ch_fcsgx,
        ch_gxdb,
        ch_ramdisk
    )

    ch_cleangenome_in  = ch_fasta
                       | join(FCSGX_RUNGX.out.fcsgx_report, by: 0)
                       | map { meta, fasta, fcsgx_report -> tuple(meta, fasta, fcsgx_report) }

    FCSGX_CLEANGENOME(
        ch_cleangenome_in
    )

    // Not cleaned
    TIARA_RAW(
        ch_fasta
    )

    // Cleaned genome
    TIARA_CLEANED(
        FCSGX_CLEANGENOME.out.cleaned
    )

    emit:
    fcs_gx_report    = FCSGX_RUNGX.out.fcsgx_report         // channel: [ val(meta), path(*.fcs_gx_report.txt) ]
    adaptor_report   = FCS_FCSADAPTOR.out.adaptor_report     // channel: [ val(meta), path(*.fcs_adaptor_report.txt) ]
    tiara_raw        = TIARA_RAW.out.classifications          // channel: [ val(meta), path(*.txt) ]
    tiara_cleaned    = TIARA_CLEANED.out.classifications      // channel: [ val(meta), path(*.txt) ]

}
