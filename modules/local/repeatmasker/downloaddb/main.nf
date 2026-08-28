process REPEATMASKER_DOWNLOADDB {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/curl:7.80.0'
        : 'biocontainers/curl:7.80.0' }"

    input:
    tuple val(meta), val(db_url)

    output:
    tuple val(meta), path("*.h5"), emit: h5
    tuple val("${task.process}"), val('<curl>'), eval('curl --version | head -1 | cut -d " " -f2'), emit: versions, topic: versions

    script:
    def filename     = db_url.tokenize('/')[-1]
    def md5_url      = db_url + ".md5"

    """
    # Download a RepeatMasker/Dfam database file and verify it against its MD5 checksum
    curl -O "${db_url}"
    curl -O "${md5_url}"

    expected_md5=\$(cat ${filename}.md5 | cut -d ' ' -f1)
    actual_md5=\$(md5sum ${filename} | cut -d ' ' -f1)

    if [ "\${actual_md5}" == "\${expected_md5}" ]; then
        gunzip ${filename}
    else
        echo "ERROR: MD5 mismatch for ${filename}"
        exit 1
    fi
    """

    stub:
    def filename   = db_url.tokenize('/')[-1]
    def out_gunzip = filename.replaceAll(/\.gz$/, '')

    """
    touch ${out_gunzip}
    """
}
