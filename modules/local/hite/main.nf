process HITE {
    tag "$meta.id"
    label 'process_medium'

    container "docker.io/kanghu/hite:3.3.3"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*_hite_results")           , emit: hite_results
    tuple val(meta), path("*_hite_results/*.tbl")  , emit: tbl
    tuple val("${task.process}"), val('python'), eval('python --version | cut -f 2 -d " "'), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('hite'), val('3.3.3'), emit: versions_hite, topic: versions
    tuple val("${task.process}"), val('repeatmasker'), eval('RepeatMasker | grep version | cut -f 3 -d " "'), emit: versions_repeatmasker, topic: versions
    tuple val("${task.process}"), val('repeatmodeler'), eval('RepeatModeler | grep /opt/conda/envs/HiTE/share/RepeatModeler/RepeatModeler | cut -f 3 -d " "'), emit: versions_repeatmodeler, topic: versions
    tuple val("${task.process}"), val('ltrpipeline'), eval('LTRPipeline -version'), emit: versions_ltrpipeline, topic: versions

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "HITE module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Run HiTE to detect and classify transposable elements in the genome assembly
    # Unzip the genome and make sure it does not have internal new line characters.
    if [ -f *.gz ]; then
      gunzip -c "$fasta" > myunzip.fa
      #myunzip.fa=\$(gunzip -c "$fasta")
      awk '/^>/ { print (NR==1 ? "" : RS) \$0; next } { printf "%s", \$0 } END { printf RS }' myunzip.fa > ${prefix}.unwrapped.fasta
    else
      awk '/^>/ { print (NR==1 ? "" : RS) \$0; next } { printf "%s", \$0 } END { printf RS }' $fasta > ${prefix}.unwrapped.fasta
    fi

    # Capture the current working directory
    mydir=`pwd`

    # Create the output directory
    mkdir -p \${mydir}/${prefix}_hite_results

    newpath=`realpath ${prefix}.unwrapped.fasta`

    cd /HiTE

    work_dir=\${TMPDIR:-/tmp}

    python main.py \\
    --genome \${newpath} \\
    --out_dir \${mydir}/${prefix}_hite_results \\
    --thread ${task.cpus} \\
    --work_dir \${work_dir} \\
    $args

    # HiTE exits 0 even when it finds no transposable elements (small/low-repeat
    # genomes), in which case it never writes HiTE.tbl. Treat that as a valid
    # zero-TE result: an empty .tbl parses as 100% non-repeat downstream.
    if [ ! -f \${mydir}/${prefix}_hite_results/HiTE.tbl ]; then
        touch \${mydir}/${prefix}_hite_results/HiTE.tbl
    fi

    mv \${mydir}/${prefix}_hite_results/HiTE.tbl \${mydir}/${prefix}_hite_results/${prefix}.tbl
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}_hite_results
    touch ${prefix}_hite_results/${prefix}.tbl
    """
}
