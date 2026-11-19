process MASH_TRIANGLE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mash:2.3--he348c14_1' :
        'quay.io/biocontainers/mash:2.3--he348c14_1' }"

    input:
    tuple val(meta), val(ids), path(fastas, stageAs: 'genomes/?/*') // ids and fastas in the same order

    output:
    tuple val(meta), path("*.tsv"), emit: distances
    tuple val(meta), path("*.msh"), emit: sketch
    tuple val("${task.process}"), val('mash'), eval("mash --version 2>&1"), emit: versions_mash, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Link each genome under its sample id, so the sketch (and the tree) is labelled by sample
    def links  = [ ids, fastas ].transpose().collect { id, fasta -> "ln -s ../${fasta} named/${id}" }.join('\n    ')
    """
    # Sketch every genome and compute all pairwise Mash distances (lower triangle)
    mkdir named
    ${links}

    cd named
    mash sketch \\
        ${args} \\
        -p ${task.cpus} \\
        -o ../${prefix} \\
        *
    cd ..

    mash triangle \\
        -p ${task.cpus} \\
        ${prefix}.msh \\
        > ${prefix}_distances.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.msh
    touch ${prefix}_distances.tsv
    """
}
