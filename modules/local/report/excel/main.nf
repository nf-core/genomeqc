process REPORT_EXCEL {
    tag "genomeqc_excel"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.9--1'
        : 'quay.io/biocontainers/python:3.9--1'}"

    input:
    path busco_tables,        stageAs: "busco/*"      // genome BUSCO batch_summary_modified.txt files
    path busco_prot_tables,   stageAs: "busco_prot/*" // protein BUSCO batch_summary_modified.txt files (optional)
    path quast_tsvs,          stageAs: "quast/*"      // per-species QUAST TSV files
    path agat_stats,          stageAs: "agat/*"       // AGAT spstatistics *.txt files
    path tidk_tsvs,           stageAs: "tidk/*"       // tidk aposteriori search TSV files
    path fcs_gx_reports,      stageAs: "fcsgx/*"      // FCS-GX *.fcs_gx_report.txt (optional)
    path fcs_adaptor_reports, stageAs: "fcsadaptor/*" // FCS-Adaptor *.fcs_adaptor_report.txt (optional)
    path tiara_reports,       stageAs: "tiara/*"      // Tiara classification *.txt (optional)
    path busco_seqs_table,    stageAs: "busco_seqs/*" // ortho_seqs.py output (optional)
    path repeatmasker_tbls,   stageAs: "repeatmasker/*" // RepeatMasker *.tbl files (optional)

    output:
    path "genomeqc_tables.xlsx", emit: excel
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //g"'), emit: versions_python, topic: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def busco_arg      = busco_tables        ? "--busco_tables busco/*"              : ""
    def busco_prot_arg = busco_prot_tables   ? "--busco_prot_tables busco_prot/*"    : ""
    def quast_arg      = quast_tsvs          ? "--quast_tsvs quast/*"                : ""
    def agat_arg       = agat_stats          ? "--agat_stats agat/*"                 : ""
    def tidk_arg       = tidk_tsvs           ? "--tidk_tsvs tidk/*"                  : ""
    def fcsgx_arg      = fcs_gx_reports      ? "--fcs_gx_reports fcsgx/*"            : ""
    def fcsadp_arg     = fcs_adaptor_reports ? "--fcs_adaptor_reports fcsadaptor/*"  : ""
    def tiara_arg      = tiara_reports       ? "--tiara_reports tiara/*"             : ""
    def busco_seqs_arg = busco_seqs_table    ? "--busco_seqs_table busco_seqs/*"     : ""
    def repeatmasker_arg = repeatmasker_tbls ? "--repeatmasker_tbls repeatmasker/*"  : ""

    """
    # Compile all per-tool QC tables into a single multi-sheet Excel workbook
    generate_excel.py \\
        ${busco_arg} \\
        ${busco_prot_arg} \\
        ${quast_arg} \\
        ${agat_arg} \\
        ${tidk_arg} \\
        ${fcsgx_arg} \\
        ${fcsadp_arg} \\
        ${tiara_arg} \\
        ${busco_seqs_arg} \\
        ${repeatmasker_arg} \\
        --output genomeqc_tables.xlsx
    """

    stub:

    """
    touch genomeqc_tables.xlsx
    """
}
