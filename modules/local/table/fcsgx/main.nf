process TABLE_FCSGX {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.9--1'
        : 'quay.io/biocontainers/python:3.9--1'}"

    input:
    tuple val(meta), path(reports)

    output:
    tuple val(meta), path("*.tsv"), emit: table
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //g"'), emit: versions_python, topic: versions


    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # Parse FCS-GX genome reports into a combined contamination-summary TSV
    fcsgx_report_2_table.py \\
        $reports \\
        -o ${prefix}.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.tsv
    """
}
