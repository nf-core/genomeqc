process REPEATMODELER_REPEATMODELER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/repeatmodeler:2.0.7--pl5321hdfd78af_0':
        'biocontainers/repeatmodeler:2.0.7--pl5321hdfd78af_0' }"

    input:
    tuple val(meta), path(db)

    output:
    tuple val(meta), path("*.fa") , emit: fasta, optional: true
    tuple val(meta), path("*.stk"), emit: stk, optional: true
    tuple val(meta), path("*.log"), emit: log, optional: true
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def db_name = file(db[0]).getBaseName()
    """
    RepeatModeler \\
        -database $db_name \\
        $args \\
        -threads $task.cpus
    
    if [ -f "${db_name}-families.fa" ]
    then
      mv ${db_name}-families.fa   ${prefix}.fa
      mv ${db_name}-families.stk  ${prefix}.stk
      mv ${db_name}-rmod.log      ${prefix}.log
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repeatmodeler: \$(RepeatModeler --version | sed 's/RepeatModeler version //')
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fa
    touch ${prefix}.stk
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repeatmodeler: \$(RepeatModeler --version | sed 's/RepeatModeler version //')
    END_VERSIONS
    """
}
