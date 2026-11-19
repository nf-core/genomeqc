process NJTREE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/77/775cbc98f9e32d9e07c11116df02f3f0f23c94c1001fa8a38825fb0183009e0c/data'
        : 'community.wave.seqera.io/library/python_pandas_r-base_bioconductor-ggtreeextra_pruned:5327fed29a6ac09f'}"

    input:
    tuple val(meta), path(distances) // mash triangle lower-triangle distance matrix

    output:
    tuple val(meta), path("*.nwk"), emit: tree
    tuple val("${task.process}"), val('r_base'), eval('Rscript -e "cat(as.character(getRversion()))"'), emit: versions_r_base, topic: versions
    tuple val("${task.process}"), val('ape'), eval('Rscript -e "cat(as.character(packageVersion(\'ape\')))"'), emit: versions_ape, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}_tree"
    """
    # Neighbour-joining tree from the distances, rooted at its midpoint
    nj_tree.R ${distances} ${prefix}.nwk
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_tree"
    """
    touch ${prefix}.nwk
    """
}
