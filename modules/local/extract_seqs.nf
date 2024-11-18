process EXTRACT_SEQS {
    tag "$meta.id"
    label 'process_low'
    //label 'process_med_memory'

    container = 'community.wave.seqera.io/library/agat:1.4.1--304a47c62ae478b4'

    input:
    tuple val (meta),  path(fasta)
    tuple val (meta),  path(gff)

    output:
    tuple val (meta), path( "${meta.id}.prot.fasta" ), emit: prot_fasta
    path "versions.yml"                              , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Extract protein from filtered GFF and genome
    agat_sp_extract_sequences.pl \\
    -g ${gff} \\
    -f ${fasta} \\
    -p -o ${prefix}.prot.fasta --clean_final_stop
    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Perl version: \$(perl --version | grep "version" | sed 's/.*(//g' | sed 's/[)].*//')
    END_VERSIONS
    """

}