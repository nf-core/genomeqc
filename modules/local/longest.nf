process LONGEST {

    tag "$meta.id"
    label 'process_single'
    label 'process_med_memory'

    container = 'biocontainers/agat:1.3.0--pl5321hdfd78af_0'
    publishDir "$params.outdir/output_data/longest" , mode: "${params.publish_dir_mode}", pattern:"*.txt"

    input:
    tuple val (meta),  path(gff)

    output:
    tuple val (meta), path( "${meta.id}.longest.gff3" ),                emit: longest_proteins
    tuple val (meta), path( "${meta.id}.stat.original.txt" ),           emit: agat_summary_original
    tuple val (meta), path( "${meta.id}.stat.long.txt" ),               emit: agat_summary_longest
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Run agat to find longest orf for each gene 
    agat_sp_keep_longest_isoform.pl -gff ${gff} -o ${prefix}.longest.gff3
    
    # Run a few summarisation scripts to report the actual genes being considered.
    agat_sp_functional_statistics.pl --gff ${gff} -o ${prefix}.stat.original.txt
    agat_sp_functional_statistics.pl --gff ${prefix}.longest.gff3 -o ${prefix}.stat.long.txt
    
    md5sum "${prefix}.longest.gff3" > "${prefix}.longest.gff3.md5"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Perl version: \$(perl --version | grep "version" | sed 's/.*(//g' | sed 's/[)].*//')
    END_VERSIONS
    """

}
