include { NCBIGENOMEDOWNLOAD                       } from '../../../modules/nf-core/ncbigenomedownload/main'
include { PIGZ_UNCOMPRESS as UNCOMPRESS_FASTA      } from '../../../modules/nf-core/pigz/uncompress/main'
include { PIGZ_UNCOMPRESS as UNCOMPRESS_GXF        } from '../../../modules/nf-core/pigz/uncompress/main'
include { BUSCO_DOWNLOAD                           } from '../../../modules/nf-core/busco/download/main'


workflow INPUT_PREPARATION {

    take:
    ch_ncbi                 // channel: [ val(meta), val(accession), val(fastq) ]
    ch_local                // channel: [ val(meta), val(fasta), val(gxf), val(fastq) ]
    val_groups               // val: NCBI genome download assembly groups (e.g. 'all')
    val_busco_lineages_path // val: path to a pre-staged BUSCO lineages directory, or null
    val_busco_lineage       // val: BUSCO lineage name (e.g. 'hymenoptera_odb10'), 'auto', or null

    main:

    ch_ncbi_meta        = ch_ncbi
                        | map { meta, _accession, _fq -> tuple( meta.id, meta ) }
    ch_ncbi_accessions  = ch_ncbi
                        | map { meta, accession, _fq -> tuple( meta.id, accession ) }
                        | collectFile { id, accession_n -> [ "${id}.txt", "${accession_n}\n" ] }
                        | map { accessions_file -> tuple( accessions_file.baseName, accessions_file ) }

    ch_ncbi_input = ch_ncbi_meta
                    | join(ch_ncbi_accessions) // join by id (meta.id and accessions_file.baseName are the same)
                    | multiMap {
                        _id, meta, accessions_file ->
                            meta      : meta
                            accession : accessions_file
                    }

    //
    // MODULE: Run ncbigenomedownload for RefSeq IDs
    //

    NCBIGENOMEDOWNLOAD (
        ch_ncbi_input.meta,
        ch_ncbi_input.accession,
        [],
        val_groups
    )

    //
    // Prepare fasta channels
    //

    // fasta. We use mix() here because when local files are present,
    // then RefSeq IDs should be missing, and viceversa
    fasta        = ch_local
                 | map { meta, fasta, _gxf, _fq -> tuple( meta, fasta) }
                 | mix ( NCBIGENOMEDOWNLOAD.out.fna )

    // Filter fasta files by extension and create channels for each file type
    gz_fasta     = fasta.filter { _meta, fasta_compressed -> fasta_compressed.name.endsWith(".gz") }
    non_gz_fasta = fasta.filter { _meta, fasta_uncompressed -> !fasta_uncompressed.name.endsWith(".gz") }

    UNCOMPRESS_FASTA ( gz_fasta )
    ch_fasta     = UNCOMPRESS_FASTA.out.file.mix(non_gz_fasta)

    //
    // Prepare gxf channels
    //

    // gxf. We use mix() here because when local files are present,
    // then RefSeq IDs should be missing, and viceversa
    gxf         = ch_local
                | map { meta, _fasta, gxf, _fq ->  tuple( meta,  gxf) }
                | mix ( NCBIGENOMEDOWNLOAD.out.gff )

    // Filter gxf files by extension and create channels for each file type
    gz_gxf      = gxf.filter { _meta, gxf_compressed -> gxf_compressed  && gxf_compressed.name.endsWith(".gz")  }
    non_gz_gxf  = gxf.filter { _meta, gxf_uncompressed -> !gxf_uncompressed || !gxf_uncompressed.name.endsWith(".gz") }

    UNCOMPRESS_GXF( gz_gxf )
    ch_gxf      = UNCOMPRESS_GXF.out.file.mix(non_gz_gxf)

    //
    // Prepare fastq channel
    //

    // FASTQ file is optional in the samplesheet.
    ch_fastq = ch_ncbi
                | map{ meta, _refseq, fq -> tuple( meta, fq ) }
                | mix( ch_local.map { meta, _fasta_local, _gxf_local, fq -> tuple( meta, fq ) } )

    //
    // Define multi-channel objects for every process/subworkflow
    //

    // Combine fasta, gxf and fastq into a single joined channel so they are in sync
    ch_input       = ch_fasta
                   | join(ch_gxf, remainder: true)    // full outer: gxf is optional
                   | join(ch_fastq, remainder: true)  // left outer: fastq is optional

    // Decontamination: filter rows where taxid is present
    ch_input_decon = ch_fasta.filter { meta, _fasta_decon -> meta.taxid }

    //
    // MODULE: Download BUSCO database
    //

    if ( val_busco_lineages_path ) {
        ch_busco_db = val_busco_lineages_path
    } else if ( val_busco_lineage != 'auto' && val_busco_lineage != null ) {
        BUSCO_DOWNLOAD ( val_busco_lineage )
        ch_busco_db = BUSCO_DOWNLOAD.out.download_dir
    } else {
        ch_busco_db = channel.value([])
    }

    emit:
    fasta       = ch_fasta        // channel: [ val(meta), path(fasta) ]
    input       = ch_input        // channel: [ val(meta), path(fasta), path(gxf), path(fastq) ]
    input_decon = ch_input_decon  // channel: [ val(meta), path(fasta) ]
    busco_db    = ch_busco_db     // channel: path(db) or val([])

}
