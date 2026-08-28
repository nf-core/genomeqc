process BUSCO_TSVTOGFF {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.9--1'
        : 'quay.io/biocontainers/python:3.9--1'}"

    input:
    tuple val(meta), path(busco_dir)

    output:
    tuple val(meta), path("*_busco.gff"), emit: gff
    tuple val(meta), path("*_busco_stats.json"), emit: stats
    //tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //g"'), emit: versions_python, topic: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    # Convert a BUSCO full_table.tsv into a GFF3 of gene locations plus a JSON stats summary
    busco_tsv_to_gff.py ${busco_dir} ${prefix}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_busco.gff
    touch ${prefix}_busco_stats.json
    """
}
