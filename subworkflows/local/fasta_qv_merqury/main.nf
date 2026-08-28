include { MERYL_UNIONSUM                           } from '../../../modules/nf-core/meryl/unionsum/main'
include { MERYL_COUNT                              } from '../../../modules/nf-core/meryl/count/main'
include { MERQURY_MERQURY                          } from '../../../modules/nf-core/merqury/merqury/main'



workflow FASTA_QV_MERQURY {

    take:
    ch_fasta
    ch_fastq
    kvalue

    main:

    // MODULE: MERYL_COUNT
    MERYL_COUNT(
        ch_fastq,
        kvalue
    )
    ch_meryl_db = MERYL_COUNT.out.meryl_db
    // MODULE: MERYL_UNIONSUM
    MERYL_UNIONSUM(
        ch_meryl_db,
        kvalue
    )
    ch_meryl_union = MERYL_UNIONSUM.out.meryl_db

    // MODULE: MERQURY_MERQURY
    ch_merqury_inputs = ch_meryl_union.join(ch_fasta)

    MERQURY_MERQURY ( ch_merqury_inputs )

    emit:
    ch_merqury_qv                           = MERQURY_MERQURY.out.assembly_qv
    ch_merqury_stats                        = MERQURY_MERQURY.out.stats
    ch_merqury_spectra_cn_fl_png            = MERQURY_MERQURY.out.spectra_cn_fl_png
    ch_merqury_spectra_asm_fl_png           = MERQURY_MERQURY.out.spectra_asm_fl_png
    ch_hapmers_blob_png                     = MERQURY_MERQURY.out.hapmers_blob_png

}
