process SHINYAPP {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/ubuntu:22.04'
        : 'nf-core/ubuntu:22.04' }"

    input:
    tuple val(meta), path(tables), path(tree)
    path(functions)
    path(app)
    val(container_engine)


    output:
    tuple val(meta), path("app/shiny_app.sh"), emit: shiny_app
    path("app"), emit: shiny_app_files
    tuple val("${task.process}"), val("cat"), eval("cat --version 2>&1 | head -n 1 | sed 's/^.*coreutils) //; s/ .*\$//'"), emit: versions_shiny, topic: versions

    script:
    def resolved_container_engine = container_engine ?: 'docker'
    // Fully qualify the registry: this is a plain `docker run`, so (unlike the
    // Nextflow-managed processes) it does not get the docker.registry='quay.io'
    // prefix and would otherwise resolve to Docker Hub, where the image is not published.
    def docker_url                = 'community.wave.seqera.io/library/python_pandas_r-base_bioconductor-ggtreeextra_pruned:5327fed29a6ac09f'

    """
    # Package the QC tables and tree into a launchable Shiny app container script
    mkdir app

    # Mount directory is executable directory (\$0)
    echo 'cd \$(dirname \$0)' > shiny_app.sh

    # Using port 8000 as is the one usually available

    echo 'CONTAINER_ID=\$($resolved_container_engine run -d -v \$(pwd):/app -w /app -p 8000:8000 $docker_url Rscript shiny_app.R)' >> shiny_app.sh
    echo 'sleep 2' >> shiny_app.sh

    # Ensure the container is stopped and port is available again once script exits
    cat <<'EOF' >> shiny_app.sh
    cleanup() {
        echo "Stopping container \$CONTAINER_ID"
        $resolved_container_engine stop "\$CONTAINER_ID" >/dev/null 2>&1
    }
    trap cleanup EXIT
    EOF

    # Make sure app can be opened in any OS (Darwin is for mac)
    cat <<'EOF' >> shiny_app.sh
    OS_TYPE=\$(uname)
    if [[ "\$OS_TYPE" == "Darwin" ]]; then
        open http://localhost:8000
    elif [[ "\$OS_TYPE" == "Linux" ]]; then
        if command -v xdg-open >/dev/null; then
            xdg-open http://localhost:8000
        else
            echo "App running at http://localhost:8000"
        fi
    elif [[ "\$OS_TYPE" == "MINGW"* || "\$OS_TYPE" == "CYGWIN"* ]]; then
        start http://localhost:8000
    else
        echo "Please open your browser at http://localhost:8000"
    fi
    # So that container doesn't close unless ctr+c
    $resolved_container_engine logs -f "\$CONTAINER_ID"
    EOF

    chmod +x shiny_app.sh

    mv *.tsv app/
    mv tree.nw app/
    mv *.R app/
    mv shiny_app.sh app/

    """

    stub:
    """
    mkdir app
    touch app/shiny_app.sh
    """
}
