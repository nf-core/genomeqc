process GFFREAD {
    tag "$meta.id"
    label 'process_single'
    container = 'ecoflowucl/gffread_python:python-3.11.9_Linux_x86_64_perl-5.36.0'

    input:
    tuple val(meta), path(fasta)
    tuple val(meta), path(gff)

    output:
    path( "${meta.id}.prot.fa" ), emit: proteins
    tuple val(meta), path("${meta.id}.prot.fa.largestIsoform.fa" ), emit: proteins_busco
    path( "${meta.id}.prot.fa.largestIsoform.fa" ), emit: longest
    path( "${meta.id}.splicedcds.fa" )
    path( "${meta.id}.splicedexons.fa" )
    path( "${meta.id}.gff_for_jvci.gff3" ), emit: gffs
    tuple val(meta), path("${meta.id}.gff_for_jvci.gff3"), emit: gffs_agat
    path( "${meta.id}_gene_alltran_list.txt" ), emit: gene_to_isoforms
    path( "${meta.id}.splicedcds.fa.nucl.longest.fa" )
    tuple val( "${meta.id}" ), path( "${fasta}" ), emit: fasta_quast
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix         = task.ext.prefix ?: "${meta.id}"
    """
    #Check if gff3 or genome file is gzipped:
    if [[ $gff == *.gz ]]; then
      zcat $gff > gff_temp
    else
      cp $gff gff_temp
    fi

    if [[ $fasta == *.gz ]]; then
      zcat $fasta > genome_temp
    else
      cp $fasta genome_temp
    fi

    #Convert Augustus gff files if found, then do gffread to print out the nucleotide files for each gene.

    head -n 1 gff_temp > tbd

    if grep -q AUGUSTUS tbd; then 
        gtf2gff.pl <gff_temp --out=${prefix}.gff_for_jvci.gff3
        #Remove lines of the GFF that have ? in the strand section, as this cannot be parsed by gffread
        awk '\$7 != "?" { print \$0 }' ${prefix}.gff_for_jvci.gff3  > ${prefix}.gff_for_jvci.noquest.gff3
        
        gffread -w ${prefix}.splicedexons.fa -g genome_temp ${prefix}.gff_for_jvci.noquest.gff3 --table @genename
        gffread -x ${prefix}.splicedcds.fa -g genome_temp ${prefix}.gff_for_jvci.noquest.gff3 --table @genename
        gffread -y ${prefix}.prot.fa -g genome_temp ${prefix}.gff_for_jvci.noquest.gff3 -F -S --table @genename

    else
        mv gff_temp ${prefix}.gff_for_jvci.gff3
        #Remove lines of the GFF that have ? in the strand section, as this cannot be parsed by gffread
        awk '\$7 != "?" { print \$0 }' ${prefix}.gff_for_jvci.gff3  > ${prefix}.gff_for_jvci.noquest.gff3
        
        gffread -w ${prefix}.splicedexons.fa -g genome_temp ${prefix}.gff_for_jvci.noquest.gff3 --table @genename
        gffread -x ${prefix}.splicedcds.fa -g genome_temp ${prefix}.gff_for_jvci.noquest.gff3 --table @genename
        gffread -y ${prefix}.prot.fa -g genome_temp ${prefix}.gff_for_jvci.noquest.gff3 -F -S --table @genename

    fi

    gff_to_genetranshash.2.pl
    prot_fasta_to_longest.pl ${prefix}.prot.fa ${prefix}_longestisoform.txt
    fasta_topIsoform.pl ${prefix}.splicedcds.fa ${prefix}_longestisoform.txt


    #This part checks if longest isoform worked, if not we will continue with all proteins into Orthofinder. Warning sent to screen.
    #Largest isoforms has content if true
    #Largest isoforms does not have content if false. Just use full protein file (could be a genome without isoforms)

    if [[ -s ${prefix}.prot.fa.largestIsoform.fa ]];then
      echo all_good
    else
      cp ${prefix}.prot.fa ${prefix}.prot.fa.largestIsoform.fa
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread version: \$(gffread --version)
        Perl version: \$(perl --version | grep "version" | sed 's/.*(//g' | sed 's/[)].*//')
    END_VERSIONS
    """
}
