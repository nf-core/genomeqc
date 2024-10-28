// TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/modules/nf-core/
//               You can also ask for help via your pull request or on the #modules channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
// TODO nf-core: Software that can be piped together SHOULD be added to separate module files
//               unless there is a run-time, storage advantage in implementing in this way
//               e.g. it's ok to have a single module for bwa to output BAM instead of SAM:
//                 bwa mem | samtools view -B -T ref.fasta
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

process FCS_CLEANGENOME {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ncbi-fcs-gx:0.5.4--h4ac6f70_1':
        'biocontainers/ncbi-fcs-gx:0.5.4--h4ac6f70_1' }"

    input:
    tuple path(assembly), path(action_report) //fcs-gx_find.out.fcs_gx_report

    output:
    path("NCBI/*.cleaned.fasta"), emit: cleaned // genome with contamination removed 
    path("NCBI/*.contam.fasta"), emit: contam // contaminated sequences in fasta format
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "FCS_FCSGX module does not support Conda. Please use Docker / Singularity / Podman instead."
        }
    def FCSGX_VERSION = '0.5.4'
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p NCBI
    gunzip ${assembly}
    gx clean-genome \\
    --input "${assembly}" \\
    --action-report "${action_report}" \\
    --output "NCBI/${meta}.cleaned.fasta" \\
    --contam-fasta-out "NCBI/${meta}.contam.fasta"
    $args
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed -e "s/Python //g")
        FCS-GX: $FCSGX_VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "FCS_FCSGX module does not support Conda. Please use Docker / Singularity / Podman instead."
        }
    def FCSGX_VERSION = '0.5.4'
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
        
    """
    mkdir -p NCBI
    touch NCBI/${meta}.cleaned.fasta
    touch NCBI/${meta}.contam.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed -e "s/Python //g")
        FCS-GX: $FCSGX_VERSION
    END_VERSIONS
    """
}