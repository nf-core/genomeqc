process REPORT_HTML {
    tag "genomeqc_report"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.9--1'
        : 'quay.io/biocontainers/python:3.9--1'}"

    input:
    path busco_tables,       stageAs: "busco/*"         // genome BUSCO batch_summary_modified.txt files (one per species)
    path busco_prot_tables,  stageAs: "busco_prot/*"    // protein BUSCO batch_summary_modified.txt files (optional)
    path agat_stats,         stageAs: "agat/*"          // AGAT sp_statistics.pl *.stats.txt files (optional)
    path quast_tsvs,         stageAs: "quast/*"         // QUAST report.tsv files (optional)
    path tidk_tsvs,          stageAs: "tidk/*"          // tidk aposteriori search TSV files
    path tidk_apriori_tsvs,  stageAs: "tidk_apriori/*"  // tidk apriori search TSV files (optional)
    path fcsgx_reports,      stageAs: "fcsgx/*"         // FCS-GX report files (optional)
    path fcsadp_reports,     stageAs: "fcsadp/*"        // FCS-Adaptor report files (optional)
    path tiara_reports,      stageAs: "tiara/*"         // Tiara classification files (optional)
    path busco_seqs_table,   stageAs: "busco_seqs/*"    // ortho_seqs.py output (optional)
    path repeatmasker_tbls,  stageAs: "repeatmasker/*"  // RepeatMasker *.tbl files (optional)

    output:
    path "genomeqc_report.html", emit: report
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //g"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def busco_arg        = busco_tables      ? "--busco_tables busco/*"             : ""
    def busco_prot_arg   = busco_prot_tables ? "--busco_prot_tables busco_prot/*"   : ""
    def agat_arg         = agat_stats        ? "--agat_stats agat/*"                : ""
    def quast_arg        = quast_tsvs        ? "--quast_tsvs quast/*"               : ""
    def tidk_tsv_arg     = tidk_tsvs         ? "--tidk_tsvs tidk/*"                 : ""
    def tidk_apriori_arg = tidk_apriori_tsvs ? "--tidk_apriori_tsvs tidk_apriori/*" : ""
    def fcsgx_arg        = fcsgx_reports     ? "--fcsgx_reports fcsgx/*"            : ""
    def fcsadp_arg       = fcsadp_reports    ? "--fcsadp_reports fcsadp/*"          : ""
    def tiara_arg        = tiara_reports     ? "--tiara_reports tiara/*"            : ""
    def busco_seqs_arg   = busco_seqs_table  ? "--busco_seqs_table busco_seqs/*"    : ""
    def repeatmasker_arg = repeatmasker_tbls ? "--repeatmasker_tbls repeatmasker/*" : ""

    """
    # Build a self-contained HTML QC report from the available per-tool tables
    generate_report.py \\
        ${busco_arg} \\
        ${busco_prot_arg} \\
        ${agat_arg} \\
        ${quast_arg} \\
        ${tidk_tsv_arg} \\
        ${tidk_apriori_arg} \\
        ${fcsgx_arg} \\
        ${fcsadp_arg} \\
        ${tiara_arg} \\
        ${busco_seqs_arg} \\
        ${repeatmasker_arg} \\
        --output genomeqc_report.html
    """

    stub:
    """
    touch genomeqc_report.html
    """
}
