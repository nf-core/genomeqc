process TREE_SUMMARY {
    tag "$meta.id"
    label 'process_single'

    container = 'ecoflowucl/genomeqc_tree:v1.0'
    publishDir "$params.outdir/output_data/longest" , mode: "${params.publish_dir_mode}", pattern:"*.txt"

    input:
    tuple val(meta), path(tree)
    path (busco)

    output:
    path( "${meta}.pdf" ),                emit: figure
    path( "versions.yml"),                emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Run summary plot
    plot_tree_summary.pl ${tree}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Perl version: \$(perl --version | grep "version" | sed 's/.*(//g' | sed 's/[)].*//')
    END_VERSIONS
    """

}
