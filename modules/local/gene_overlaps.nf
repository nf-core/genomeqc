process GENE_OVERLAPS {
    tag "$meta.id"
    label 'process_single'
    container = 'ecoflowucl/gene_overlap:v1.0'

    input:
    tuple val(meta), path(gff)

    output:
    path( "Summary.tsv" ), emit: overlap_summary

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix         = task.ext.prefix ?: "${meta.id}"
    """
    #Run overlap R script

    gene_overlaps.R $gff Summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(as.character(getRversion()))")
    END_VERSIONS
    """
}
