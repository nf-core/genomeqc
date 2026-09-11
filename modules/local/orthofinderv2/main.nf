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
    tuple val(meta), path("$results_dir")                                          , emit: orthofinder
    path("$results_dir/Phylogenetic_Hierarchical_Orthogroups/N0.tsv")              , emit: orthologues
    path("$results_dir/Orthogroups/Orthogroups.tsv")                               , emit: orthogroups
    path("$results_dir/Species_Tree/SpeciesTree_rooted_node_labels.txt")           , emit: speciestree
    tuple val("${task.process}"), val('orthofinder'), eval("orthofinder -h | sed -n 's/.*version \\(.*\\) Copy.*/\\1/p'"), emit: versions_orthofinder, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    results_dir = "input/OrthoFinder/Results_${prefix}"

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
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    results_dir = "input/OrthoFinder/Results_${prefix}"

    """
    echo $args

    mkdir -p    $results_dir/Comparative_Genomics_Statistics
    mkdir       $results_dir/Gene_Duplication_Events
    mkdir       $results_dir/Gene_Trees
    mkdir       $results_dir/Orthogroup_Sequences
    mkdir       $results_dir/Orthogroups
    mkdir       $results_dir/Orthologues
    mkdir       $results_dir/Phylogenetic_Hierarchical_Orthogroups
    mkdir       $results_dir/Phylogenetically_Misplaced_Genes
    mkdir       $results_dir/Putative_Xenologs
    mkdir       $results_dir/Resolved_Gene_Trees
    mkdir       $results_dir/Single_Copy_Orthologue_Sequences
    mkdir       $results_dir/Species_Tree
    mkdir       $results_dir/WorkingDirectory
    touch       $results_dir/Log.txt
    touch       $results_dir/Orthogroups/Orthogroups.tsv
    touch       $results_dir/Species_Tree/SpeciesTree_rooted_node_labels.txt
    touch       $results_dir/Phylogenetic_Hierarchical_Orthogroups/N0.tsv
    """
}
