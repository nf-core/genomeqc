process TABLE_TE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.9--1'
        : 'quay.io/biocontainers/python:3.9--1'}"

    input:
    tuple val(meta), path(tables)

    output:
    tuple val(meta), path("*.tsv"), emit: table
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //g"'), emit: versions_python, topic: versions


    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # Parse RepeatMasker .tbl files into a combined TE-composition summary TSV
    te_tbl_2_table.py \\
        $tables \\
        -o ${prefix}.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.tsv
    """
}
