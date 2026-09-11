process ORTHOFINDERV2 {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/orthofinder:2.5.5--hdfd78af_2'
        : 'biocontainers/orthofinder:2.5.5--hdfd78af_2' }"

    input:
    tuple val(meta), path(fastas, stageAs: 'input/')

    output:
    tuple val(meta), path("$prefix")                                          , emit: orthofinder
    path("$prefix/Phylogenetic_Hierarchical_Orthogroups/N0.tsv")              , emit: orthologues
    path("$prefix/Species_Tree/SpeciesTree_rooted_node_labels.txt")           , emit: speciestree
    tuple val("${task.process}"), val('orthofinder'), eval("orthofinder -h | sed -n 's/.*version \\(.*\\) Copy.*/\\1/p'"), emit: versions_orthofinder, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"

    """
    # Infer orthogroups and a species tree across the input proteomes with OrthoFinder2
    mkdir temp_pickle

    orthofinder \\
        $args \\
        -t $task.cpus \\
        -a ${[task.cpus, 4].min()} \\
        -p temp_pickle \\
        -f input \\
        -n $prefix

    mv \\
    input/OrthoFinder/Results_$prefix \\
    $prefix
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo $args

    mkdir -p    $prefix/Comparative_Genomics_Statistics
    mkdir       $prefix/Gene_Duplication_Events
    mkdir       $prefix/Gene_Trees
    mkdir       $prefix/Orthogroup_Sequences
    mkdir       $prefix/Orthogroups
    mkdir       $prefix/Orthologues
    mkdir       $prefix/Phylogenetic_Hierarchical_Orthogroups
    mkdir       $prefix/Phylogenetically_Misplaced_Genes
    mkdir       $prefix/Putative_Xenologs
    mkdir       $prefix/Resolved_Gene_Trees
    mkdir       $prefix/Single_Copy_Orthologue_Sequences
    mkdir       $prefix/Species_Tree
    mkdir       $prefix/WorkingDirectory
    touch       $prefix/Log.txt
    touch       $prefix/Orthogroups/Orthogroups.tsv
    touch       $prefix/Species_Tree/SpeciesTree_rooted_node_labels.txt
    touch       $prefix/Phylogenetic_Hierarchical_Orthogroups/N0.tsv
    """
}
