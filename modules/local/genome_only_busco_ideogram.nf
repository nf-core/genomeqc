process GENOME_ONLY_BUSCO_IDEOGRAM {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::r-rideogram=0.2.2 bioconda::r-svglite=2.1.1"
    container "community.wave.seqera.io/library/seqkit_r-dplyr_r-optparse_r-readr_pruned:92b750716d244919"

    input:
    tuple val(meta), path(genome), path(busco_full_table)

    output:
    tuple val(meta), path("*.svg"), emit: svg
    tuple val(meta), path("*.png"), emit: png
    path "versions.yml"           , emit: versions

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    grep -v "#" ${busco_full_table} | cut -f 2,3,4,5  | grep -v "Missing" > ${prefix}_busco_coordinates.txt

    # Get chromosome lengths:
    seqkit fx2tab -i -n -l ${genome} > ${prefix}_for_karyotype.txt

    # Call script for table wrangling
    plot_markers1.R ${prefix}_for_karyotype.txt ${prefix}

    # Call script for plotting
    plot_busco_ideogram.R \\
        --busco_output ${prefix}_busco_coordinates.txt \\
        --karyotype ${prefix}_karyotype.txt \\
        --prefix ${prefix} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(as.character(getRversion()))")
        r-rideogram: \$(Rscript -e "cat(as.character(packageVersion('RIdeogram')))")
    END_VERSIONS
    """
}
