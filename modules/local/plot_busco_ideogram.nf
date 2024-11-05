process PLOT_BUSCO_IDEOGRAM {
    tag "${genusspeci}_${lineage}"
    label 'process_low'

    conda "bioconda::r-rideogram=0.2.2 bioconda::r-svglite=2.1.1"
    container "community.wave.seqera.io/library/r-optparse_r-rideogram:14e26839d69da37f"

    input:
    tuple val(genusspeci), val(lineage), path(busco_full_table), path(gff_file)

    output:
    tuple val(genusspeci), val(lineage), path("*.svg"), emit: svg
    tuple val(genusspeci), val(lineage), path("*.png"), emit: png
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = "${genusspeci}_${lineage}"
    """
    plot_busco_ideogram.R \\
        --busco_full_table ${busco_full_table} \\
        --gff_file ${gff_file} \\
        --lineage ${lineage} \\
        --prefix ${prefix} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(as.character(getRversion()))")
        r-rideogram: \$(Rscript -e "cat(as.character(packageVersion('RIdeogram')))")
    END_VERSIONS
    """
}