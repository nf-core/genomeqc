process BUSCO_SEQS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/python_pip_pandas:2fd05a70c67560f2"

    input:
    tuple val(meta), path(tables)

    output:
    tuple val(meta), path("*.tsv"), emit: table
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix         = task.ext.prefix ?: "${meta.id}"
    """
    # Get chromosome lengths:
    ortho_seqs.py \\
    -i $tables \\
    $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
