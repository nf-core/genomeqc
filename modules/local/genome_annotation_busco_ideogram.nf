process GENOME_ANNOTATION_BUSCO_IDEOGRAM {
    tag "${genusspeci}_${lineage}"
    label 'process_single'

    conda "bioconda::r-rideogram=0.2.2 bioconda::r-svglite=2.1.1"
    container "community.wave.seqera.io/library/seqkit_r-dplyr_r-optparse_r-readr_pruned:92b750716d244919"

    input:
    tuple val(genusspeci), val(lineage), path(busco_full_table), path(genome), path(gff)

    output:
    tuple val(genusspeci), val(lineage), path("*.svg"), emit: svg
    tuple val(genusspeci), val(lineage), path("*.png"), emit: png
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = "${genusspeci}_${lineage}"
    """
    # Get chromosome lengths:
    seqkit fx2tab -i -n -l ${genome} > ${prefix}_for_karyotype.txt

    # Call the R script for plotting
    plot_busco_ideogram_kayrotyope_setup.R ${prefix}_for_karyotype.txt ${prefix}

    #Extract the information from gff and busco report to get locations of busco genes:
    busco_create_table_for_plot.R ${busco_full_table}  ${gff} busco_data_to_plot.tsv

    plot_busco_ideogram.R \\
        --busco_output busco_data_to_plot.tsv \\
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
