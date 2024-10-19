process TREE_SUMMARY {

    tag "$meta1.id"
    label 'process_single'

    container = 'biocontainers/agat:1.3.0--pl5321hdfd78af_0'
    publishDir "$params.outdir/output_data/longest" , mode: "${params.publish_dir_mode}", pattern:"*.txt"

    input:
    tuple val(meta1), path(tree)
    path (busco)

    output:
    path( "${meta1}.pdf" ),                emit: figure
    path( "versions.yml"    ),                emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta1.id}"
    """
    # Run summary plot
    plot_tree_summary.pl ${tree}

    """

}
