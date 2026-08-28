process ORTHOSEQCOUNT {
    tag "ortho_seq_count"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b3/b39f468dbf576d7e8e3e2913cb1deaaf172f194bf7c76ce80a2b8940a7c492a4/data'
        : 'community.wave.seqera.io/library/python_pip_pandas:2fd05a70c67560f2'}"

    input:
    path orthogroups_tsv
    path gff_files

    output:
    path "species_ortho_seq_count.tsv", emit: species_summary
    path "pairwise_ortho_seq_count.tsv", emit: pairwise_summary
    path "debug_gene_mapping.txt", emit: debug_info
    tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //g"'), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('pandas'), eval('python -c "import pandas as pd; print(pd.__version__)"'), emit: versions_pandas, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Map orthogroup genes onto their sequences (chromosomes/scaffolds/contigs) via the
    # per-species GFFs, then tally how often sequence pairs co-occur across species.
    ortho_seq_count.py \\
        --orthogroups ${orthogroups_tsv} \\
        --gff ${gff_files} \\
        --pairwise-out pairwise_ortho_seq_count.tsv \\
        --species-out species_ortho_seq_count.tsv \\
        --debug-out debug_gene_mapping.txt
    """

    stub:

    """
    touch species_ortho_seq_count.tsv
    touch pairwise_ortho_seq_count.tsv
    touch debug_gene_mapping.txt
    """
}
