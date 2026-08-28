process TREESUMMARY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/77/775cbc98f9e32d9e07c11116df02f3f0f23c94c1001fa8a38825fb0183009e0c/data'
        : 'community.wave.seqera.io/library/python_pandas_r-base_bioconductor-ggtreeextra_pruned:5327fed29a6ac09f'}"

    input:
    tuple val(meta), path(tree)
    tuple val(meta1), path(geno_busco)
    tuple val(meta2), path(prot_busco)
    tuple val(meta3), path(te_table)
    tuple val(meta4), path(fcs_table)
    path  multiqc_files

    output:
    path( "*.pdf"          ),                 emit: figure
    path( "*.svg"          ),                 emit: figure_svg
    tuple val(meta), path("*.tsv"),           emit: tables
    tuple val(meta), path("tree.nw"),         emit: tree
    tuple val("${task.process}"), val('r_base'), eval('Rscript -e "cat(as.character(getRversion()))"'), emit: versions_r_base, topic: versions
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //g"'), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('ggtree'), eval('Rscript -e "cat(as.character(packageVersion(\'ggtree\')))"'), emit: versions_ggtree, topic: versions
    tuple val("${task.process}"), val('ggplot2'), eval('Rscript -e "cat(as.character(packageVersion(\'ggplot2\')))"'), emit: versions_ggplot2, topic: versions
    tuple val("${task.process}"), val('patchwork'), eval('Rscript -e "cat(as.character(packageVersion(\'patchwork\')))"'), emit: versions_patchwork, topic: versions
    tuple val("${task.process}"), val('argparse'), eval('Rscript -e "cat(as.character(packageVersion(\'argparse\')))"'), emit: versions_argparse, topic: versions
    tuple val("${task.process}"), val('dplyr'), eval('Rscript -e "cat(as.character(packageVersion(\'dplyr\')))"'), emit: versions_dplyr, topic: versions
    tuple val("${task.process}"), val('tidyr'), eval('Rscript -e "cat(as.character(packageVersion(\'tidyr\')))"'), emit: versions_tidyr, topic: versions
    tuple val("${task.process}"), val('scatterpie'), eval('Rscript -e "cat(as.character(packageVersion(\'scatterpie\')))"'), emit: versions_scatterpie, topic: versions
    tuple val("${task.process}"), val('scales'), eval('Rscript -e "cat(as.character(packageVersion(\'scales\')))"'), emit: versions_scales, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args     ?: ''
    def counts_command = meta.mode == 'genome_anno' ? "gene_overlaps_table.py *.counts.tsv gene_stats.tsv --include-same-strand --include-opposite-strand" : "touch gene_stats.tsv" // Genome only option needs a gene_stats file, even if it's empty. Should check for a more elegant solution
    def ortho_file = file("species_ortho_seq_count.tsv") ? "--ortho_file species_ortho_seq_count.tsv" : ''
    def geno_busco_combined = geno_busco ? '''{ head -qn 1 *-genome-busco.batch_summary_modified.txt | head -n 1; tail -q -n 1 *-genome-busco.batch_summary_modified.txt | sed -E 's/\t+/\t/g' | sed 's/\t$//g'; } > Busco_combined_geno.tsv''' : ''
    def prot_busco_combined = prot_busco ? '''{ head -qn 1 *-proteins-busco.batch_summary_modified.txt | head -n 1; tail -q -n 1 *-proteins-busco.batch_summary_modified.txt | sed -E 's/\t+/\t/g' | sed 's/\t$//g'; } > Busco_combined_prot.tsv''' : ''
    def geno_busco_file = geno_busco ? '--busco_geno Busco_combined_geno.tsv' : ''
    def prot_busco_file = prot_busco ? '--busco_prot Busco_combined_prot.tsv' : ''
    def te_table_file = te_table ? "--te_file $te_table" : ''
    def fcs_table_file = fcs_table ? "--fcs_file $fcs_table" : ''


    """
    # Combine per-species QC tables and render the phylogenetic tree summary plot
    #Remove unwanted extensions in the tree file
    sed \'s/\\.prot\\.fa\\.largestIsoform//g\' ${tree}/Species_Tree/SpeciesTree_rooted_node_labels.txt > tree.nw

    echo $ortho_file

    # Combine GENE OVERLAPS outputs
    ${counts_command}

    # Combine the BUSCO outputs and remove empty tabs
    ${geno_busco_combined}
    ${prot_busco_combined}

    # Combine QUAST ouput
    quast_2_table.py *quast.tsv -o Quast_to_plot.tsv -col N50,N90,"Total length","GC (%)","# contigs" -plot_types bar,bar,bar,bar,bar

    rm -f *.batch_summary_modified.txt
    rm -f *.counts.tsv
    rm -f *.quast.tsv

    # Run summary plot
    plot_tree_summary.R \\
    tree.nw \\
    Quast_to_plot.tsv \\
    gene_stats.tsv \\
    n_seqs_above_x_buscos.tsv \\
    $geno_busco_file \\
    $prot_busco_file \\
    $ortho_file \\
    $te_table_file \\
    $fcs_table_file \\
    $args

    # Make sure input TSV files are captured as outputs by copying them
    # This "touches" them so Nextflow sees them as process outputs
    # This is necessary because tables are inputs to the SHINYAPP process
    if [ -f "n_seqs_above_x_buscos.tsv" ]; then
        cp n_seqs_above_x_buscos.tsv n_seqs_above_x_buscos_output.tsv
    fi

    if [ -f "species_ortho_seq_count.tsv" ]; then
        cp species_ortho_seq_count.tsv species_ortho_seq_count_output.tsv
    fi

    if [ -n "${te_table}" ] && [ -f "${te_table}" ]; then
        cp "${te_table}" te_table_output.tsv
    fi

    if [ -n "${fcs_table}" ] && [ -f "${fcs_table}" ]; then
        cp "${fcs_table}" fcs_table_output.tsv
    fi

    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo $args

    touch tree_plot.pdf
    touch tree_plot.svg
    touch tree.nw
    touch ${prefix}.tsv
    touch te_table_output.tsv
    touch fcs_table_output.tsv
    touch n_seqs_above_x_buscos_output.tsv
    touch species_ortho_seq_count_output.tsv
    """
}
