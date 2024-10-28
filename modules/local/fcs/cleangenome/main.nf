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

process CLEANGENOME {
    tag "$meta.id"
    label 'process_low'

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ncbi-fcs-gx:0.5.4--h4ac6f70_1':
        'biocontainers/ncbi-fcs-gx:0.5.4--h4ac6f70_1' }"

    input:
    // TODO nf-core: Where applicable all sample-specific information e.g. "id", "single_end", "read_group"
    //               MUST be provided as an input via a Groovy Map called "meta".
    //               This information may not be required in some instances e.g. indexing reference genome files:
    //               https://github.com/nf-core/modules/blob/master/modules/nf-core/bwa/index/main.nf
    // TODO nf-core: Where applicable please provide/convert compressed files as input/output
    //               e.g. "*.fastq.gz" and NOT "*.fastq", "*.bam" and NOT "*.sam" etc.
    tuple path(assembly), path(action_report) //fcs-gx_find.out.fcs_gx_report

    output:
    path("NCBI/*.cleaned.fasta"), emit: cleaned // genome with contamination removed 
    path("NCBI/*.contam.fasta"), emit: contam // contaminated sequences in fasta format
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: Where possible, a command MUST be provided to obtain the version number of the software e.g. 1.10
    //               If the software is unable to output a version number on the command-line then it can be manually specified
    //               e.g. https://github.com/nf-core/modules/blob/master/modules/nf-core/homer/annotatepeaks/main.nf
    //               Each software used MUST provide the software name and version number in the YAML version file (versions.yml)
    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "task.ext.args" directive
    // TODO nf-core: If the tool supports multi-threading then you MUST provide the appropriate parameter
    //               using the Nextflow "task" variable e.g. "--threads $task.cpus"
    // TODO nf-core: Please replace the example samtools command below with your module's command
    // TODO nf-core: Please indent the command appropriately (4 spaces!!) to help with readability ;)
    """
    export taxid=\$(cat "${params.results}/taxon.txt" | grep -w ${og_num} | awk -F'\\t' '{print \$4}')
    mkdir -p NCBI
    gunzip ${assembly}
    gx clean-genome \\
    --input "${sample_id}.v129mh.fasta" \\
    --action-report "${action_report}" \\
    --output "NCBI/${meta}.cleaned.fasta" \\
    --contam-fasta-out "NCBI/${meta}.contam.fasta"
                $args
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cleangenome: \$(fcs-gx --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    //               Have a look at the following examples:
    //               Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    //               Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cleangenome: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}




 input:
       
            tuple path(assembly), path(action_report) //fcs-gx_find.out.fcs_gx_report
        
        output:
            tuple val(og_num), path("NCBI/${sample_id}.v129mh.rc.fasta"), emit: cleaned
            path("NCBI/${og_num}.contam.fasta")                         , emit: contam

        script:

            // Exit if running this module with -profile conda / -profile mamba
            if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
                error "FCS_FCSGX module does not support Conda. Please use Docker / Singularity / Podman instead."
            }
            def args = task.ext.args ?: ''
            def prefix = task.ext.prefix ?: "${sample_id}"
            def FCSGX_VERSION = '0.5.4'
        
            """
            export taxid=\$(cat "${params.results}/taxon.txt" | grep -w ${og_num} | awk -F'\\t' '{print \$4}')
            mkdir -p NCBI
            gunzip ${assembly}
            gx clean-genome \\
                --input "${sample_id}.v129mh.fasta" \\
                --action-report "${action_report}" \\
                --output "NCBI/${sample_id}.v129mh.rc.fasta" \\
                --contam-fasta-out "NCBI/${og_num}.contam.fasta"
                $args
        
            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                FCS-GX: $FCSGX_VERSION
            END_VERSIONS
            """
        
            stub:
            // Exit if running this module with -profile conda / -profile mamba
            if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
                error "FCS_FCSGX module does not support Conda. Please use Docker / Singularity / Podman instead."
            }
            def prefix = task.ext.prefix ?: "${sample_id}"
            def FCSGX_VERSION = '0.5.4'
        
            """
            mkdir -p NCBI
            touch NCBI/${sample_id}.v129mh.rc.fasta
            touch NCBI/${og_num}.contam.fasta
        
            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                FCS-GX: $FCSGX_VERSION
            END_VERSIONS
            """