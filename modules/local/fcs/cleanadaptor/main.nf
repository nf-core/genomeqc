process CLEANADAPTOR {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ncbi-fcs-gx:0.5.4--h4ac6f70_1':
        'biocontainers/ncbi-fcs-gx:0.5.4--h4ac6f70_1' }"

    input:
    path(assembly), path(adaptor_report)

    output:
    tuple val(og_num), val(sample_id), path("${sample_id}.rmadapt.fasta")   , emit: cleaned
    path("${sample_id}.adaptor-contam.fasta")                               , emit: contam

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def FCSGX_VERSION = '0.5.4'

    """
    python3 --version
    #sed -i '1i ##[["FCS genome report",2,' "$adaptor_report"
    gx clean-genome \\
    --input ${filtered_fasta} \\
    --action-report "${adaptor_report}" \\
    --output "${prefix}.rmadapt.fasta" \\
    --contam-fasta-out "${prefix}.adaptor-contam.fasta"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        FCS-GX: $FCSGX_VERSION
    END_VERSIONS
    """      

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def FCSGX_VERSION = '0.5.4'

    """
    touch "${prefix}.rmadapt.fasta"
    touch "${prefix}.adaptor-contam.fasta"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        FCS-GX: $FCSGX_VERSION
    END_VERSIONS
    """
}
