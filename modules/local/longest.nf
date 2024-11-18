process LONGEST {
    tag "$meta.id"
    label 'process_medium'
    //label 'process_med_memory'

    container = 'community.wave.seqera.io/library/agat:1.4.1--304a47c62ae478b4'

    input:
    tuple val (meta),  path(gff)

    output:
    tuple val (meta), path( "${meta.id}.longest.gff3" ),                emit: longest_proteins
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Run agat to find longest orf for each gene 
    agat_sp_keep_longest_isoform.pl -gff ${gff} -o ${prefix}.longest.gff3
    
    md5sum "${prefix}.longest.gff3" > "${prefix}.longest.gff3.md5"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Perl version: \$(perl --version | grep "version" | sed 's/.*(//g' | sed 's/[)].*//')
    END_VERSIONS
    """

}
