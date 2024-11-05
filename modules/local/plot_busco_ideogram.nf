process PLOT_BUSCO_IDEOGRAM {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::r-rideogram=0.2.2 bioconda::r-svglite=2.1.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-e11b2d0c9a6d8c1b9c1c9a6d8c1b9c1c:4' :
        'biocontainers/mulled-v2-e11b2d0c9a6d8c1b9c1c9a6d8c1b9c1c:4' }"

    input:
    tuple val(meta), path(busco_full_table), path(gff_file)

    output:
    tuple val(meta), path("*.svg"), emit: svg
    tuple val(meta), path("*.png"), emit: png
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    plot_busco_ideogram.R \\
        --busco_full_table ${busco_full_table} \\
        --gff_file ${gff_file} \\
        --lineage ${meta.lineage} \\
        --prefix ${prefix} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(as.character(getRversion()))")
        r-rideogram: \$(Rscript -e "cat(as.character(packageVersion('RIdeogram')))")
    END_VERSIONS
    """
}