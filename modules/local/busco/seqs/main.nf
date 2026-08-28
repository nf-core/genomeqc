process BUSCO_SEQS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b3/b39f468dbf576d7e8e3e2913cb1deaaf172f194bf7c76ce80a2b8940a7c492a4/data'
        : 'community.wave.seqera.io/library/python_pip_pandas:2fd05a70c67560f2'}"

    input:
    tuple val(meta), path(tables)

    output:
    tuple val(meta), path("*.tsv"), emit: table
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //g"'), emit: versions_python, topic: versions

    script:
    def args = task.ext.args ?: ''

    """
    # Count sequences with Complete_BUSCOs above the threshold
    ortho_seqs.py \\
    -i $tables \\
    $args

    """

    stub:
    """
    touch n_seqs_above_x_buscos.tsv
    """
}
