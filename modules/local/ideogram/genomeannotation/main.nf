process IDEOGRAM_GENOMEANNOTATION {
    tag "${meta.id}_${meta.lineage}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d4/d4e58ef54f830d8efd697670b69953ef30b5165b812c1e0e80943e075af4b106/data'
        : 'community.wave.seqera.io/library/r-base_seqkit_r-rideogram_r-optparse_pruned:381d161bd6bee823' }"

    input:
    tuple val(meta), path(busco_full_table), path(genome), path(gff)

    output:
    tuple val(meta), path("*.svg"), emit: svg
    tuple val(meta), path("*.png"), emit: png
    tuple val(meta), path("*.csv"), emit: busco_mappings
    tuple val("${task.process}"), val('r_base'), eval('Rscript -e "cat(as.character(getRversion()))"'), emit: versions_r_base, topic: versions
    tuple val("${task.process}"), val('r_rideogram'), eval('Rscript -e "cat(as.character(packageVersion(\'RIdeogram\')))"'), emit: versions_rideogram, topic: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = "${meta.id}"

    """
    # Plot a chromosome ideogram of annotated BUSCO gene locations from the GFF
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

    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = "${meta.id}"

    """
    echo $args

    touch ${prefix}.svg
    touch ${prefix}.png
    touch ${prefix}.csv
    """
}
