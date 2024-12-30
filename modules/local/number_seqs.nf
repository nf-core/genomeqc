process NUMBER_SEQS {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::coreutils=8.31"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h36e9172_9' :
        'biocontainers/gnu-wget:1.18--h36e9172_9' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path( "${meta.id}.n_seqs.tsv" ), emit: n_seqs

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # Add header
    echo -e "species\tchromosome_count\tsequence_count" > ${prefix}.n_seqs.tsv

    # Use grep to cunt the number of chromosomes and sequences
    echo -e "${prefix}\t\$(grep -i '^>.*\\(chromosome\\|chr\\)' $fasta | wc -l)\t\$(grep -c '^>' $fasta)" >> ${prefix}.n_seqs.tsv
    """

}