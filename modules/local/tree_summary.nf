process TREE_SUMMARY {
    tag "$meta.id"
    label 'process_single'

    container = 'ecoflowucl/genomeqc_tree:v1.3'
    publishDir "$params.outdir/tree_plots" , mode: "${params.publish_dir_mode}", pattern:"Phyloplot_*.pdf"

    input:
    tuple val(meta), path(tree)
    path  multiqc_files

    output:
    path( "P*.pdf"          ),                emit: figure
    path( "versions.yml"    ),                emit: versions


    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    #Remove unwanted extensions in the tree file
    sed \'s/\\.prot\\.fa\\.largestIsoform//g\' ${tree}/Species_Tree/SpeciesTree_rooted_node_labels.txt > tree.nw
    
    # Combine the BUSCO outputs and remove empty tabs
    head -qn 1 *.txt | head -n 1                               > Busco_combined
    tail -q -n 1 *.txt | sed -E 's/\t+/\t/g' | sed 's/\t\$//g' >> Busco_combined
    # cut -f 1,3,4,5,6,7 Busco_combined >> Busco_combined_cut
    # sed -i \'s/\\.fasta//g\' Busco_combined_cut

    # busco_2_table.py Busco_combined_cut Busco_to_plot.tsv

    # Combine QUAST ouput
    quast_2_table.py *quast.tsv -o Quast_to_plot.tsv -col N50,N90 -plot_types bar,bar
    
    sudo apt-get install libmagick++-dev

    # Run summary plot
    tree_summary.R tree.nw Busco_combined Quast_to_plot.tsv
    #plot_tree_summary.R  tree.nw Quast_to_plot_final.tsv --tree_size 0.6

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Perl version: \$(perl --version | grep "version" | sed 's/.*(//g' | sed 's/[)].*//')
        Python version: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """

}
