process FAMDBPYEMBL {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/repeatmasker:4.2.2--pl5321hdfd78af_0'
        : 'biocontainers/repeatmasker:4.2.2--pl5321hdfd78af_0' }"

    input:
    tuple val(meta), path(h5_dir)
    val(lineage)

    output:
    tuple val(meta), path("*.fasta"), emit: famdb_lib
    tuple val("${task.process}"), val('repeatmasker'), eval('RepeatMasker | grep version | cut -f 3 -d " "'), emit: versions_repeatmasker, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def term = lineage ?: 'root'
    def famdb_script = workflow.containerEngine ?
        "/usr/local/share/RepeatMasker/famdb.py" :
        "\$CONDA_PREFIX/share/RepeatMasker/famdb.py"

    """
    # Export the Dfam repeat family library for a lineage and convert it to FASTA
    python ${famdb_script} \\
        -i ./ \\
        families --descendants $term -f embl \\
        $args \\
        | famdb_embl_to_fasta.py \\
        > ${term}.fasta
    """

    stub:
    def term = lineage ?: 'root'

    """
    touch ${term}.fasta
    """
}
