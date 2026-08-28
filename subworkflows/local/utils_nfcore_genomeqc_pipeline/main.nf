//
// Subworkflow with functionality specific to the nf-core/genomeqc pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { paramsHelp                } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //

    def before_text = ""
    def after_text = ""
    before_text = """
-\033[2m----------------------------------------------------\033[0m-
                                        \033[0;32m,--.\033[0;30m/\033[0;32m,-.\033[0m
\033[0;34m        ___     __   __   __   ___     \033[0;32m/,-._.--~\'\033[0m
\033[0;34m  |\\ | |__  __ /  ` /  \\ |__) |__         \033[0;33m}  {\033[0m
\033[0;34m  | \\| |       \\__, \\__/ |  \\ |___     \033[0;32m\\`-._,-`-,\033[0m
                                        \033[0;32m`._,._,\'\033[0m
\033[0;35m  nf-core/genomeqc ${workflow.manifest.version}\033[0m
-\033[2m----------------------------------------------------\033[0m-
"""
    after_text = """${workflow.manifest.doi ? "\n* The pipeline\n" : ""}${workflow.manifest.doi.tokenize(",").collect { doi -> "    https://doi.org/${doi.trim().replace('https://doi.org/','')}"}.join("\n")}${workflow.manifest.doi ? "\n" : ""}
* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
    https://github.com/nf-core/genomeqc/blob/master/CITATIONS.md
"""
    if (monochrome_logs) {
        before_text = before_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    //
    // Create channel from input file provided through params.input
    //

    channel
        .fromList(samplesheetToList(input, "${projectDir}/assets/schema_input.json"))
        .set { ch_samplesheet }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)

    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    def hite = params.te == 'hite'
    def profile_conda = workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    if (hite && profile_conda) {
        error "The HiTE module does not support Conda. Please use Docker / Singularity / Podman instead or change the TE detection method to 'repeatmasker'"
    }
    // Add ways of validating input parameters, e.g. for groups in ncbigenomedownload
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def ( meta, ncbi, fasta, gff, fastq ) = input
    // If input are ncbi accessions
    if ( meta && ncbi ) {
        if ( ncbi.startsWith('GCF') ) { // For refseq IDs
            return [ [id:meta.id, ncbi:'refseq',taxid:meta.taxid], ncbi, fastq ]
        } else if ( ncbi.startsWith('GCA') ) { // For genbank IDs
            return [ [id:meta.id, ncbi:'genbank',taxid:meta.taxid], ncbi, fastq ]
        } else {
            error('Incorrect ncbi ID. Please make sure ncbi IDs start with "GCA" for GenBank or "GCF" for RefSeq')
        }
    // If input are local files
    } else if ( meta && fasta ) { // At least fasta file is necessary if local files (genome only mode is the minimum run)
        return [ meta, fasta, gff, fastq ]
    } else {
        error("You are running genomeqc on default mode. Please check input samplesheet -> Incorrect samplesheet format")
    }
}

// MultiMap channel to split the input samplesheet into multiple channels for different processes.
// This is necessary since some processes only require a subset of the input files
// (e.g. fasta and gff for annotation-based QC, while fasta only for genome-only QC)
def multimapChannel(input) {
   def multi_ch = input
            .multiMap {
                meta, fasta, gxf, fq ->
                    fasta : fasta ? tuple( meta, file(fasta) ) : null
                    gxf   : gxf   ? tuple( meta, file(gxf) )   : null
                    fq    : fq    ? tuple( meta, file(fq) )    : null // Only this one is necessary
            }
    return multi_ch
}

//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // Tools that always run are listed unconditionally; optional tools are gated
    // on the same params that control whether they execute in the workflow.
    def citation_text = [
            "Tools used in the workflow included:",
            "FastQC (Andrews 2010),",
            "SeqKit (Shen et al. 2016),",
            "ncbi-genome-download (Blin),",
            params.genome_only ? "" : "AGAT (Dainat),",
            params.genome_only ? "" : "GffRead (Pertea and Pertea 2020),",
            "BUSCO (Manni et al. 2021),",
            "QUAST (Gurevich et al. 2013),",
            "Merqury (Rhie et al. 2020),",
            !params.skip_tidk ? "tidk (Brown et al. 2025)," : "",
            (params.gxdb || params.gxdb_manifest) ? "FCS-GX and FCS-adaptor (Astashyn et al. 2024)," : "",
            (params.gxdb || params.gxdb_manifest) ? "Tiara (Karlicki et al. 2022)," : "",
            params.ortho_version == 'v2' ? "OrthoFinder (Emms and Kelly 2019)," : "OrthoFinder (Emms et al. 2026),",
            "RIdeogram (Hao et al. 2020),",
            "ggtree (Yu et al. 2017),",
            "ape (Paradis and Schliep 2019),",
            "ggplot2 (Wickham 2016),",
            "GenomicRanges (Lawrence et al. 2013),",
            "the tidyverse (Wickham et al. 2019),",
            "pandas (The pandas development team),",
            "and MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // Entries are gated on the same params as toolCitationText() so the
    // bibliography matches the in-text citation list.
    def reference_text = [
            "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
            "<li>Shen W, Le S, Li Y, Hu F. (2016) SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLoS One, 11(10), e0163962. doi: 10.1371/journal.pone.0163962</li>",
            "<li>Blin K. ncbi-genome-download: Scripts to download genomes and metadata from the NCBI FTP servers. URL: https://github.com/kblin/ncbi-genome-download</li>",
            params.genome_only ? "" : "<li>Dainat J. AGAT: Another Gff Analysis Toolkit to handle annotations in any GTF/GFF format. Zenodo. doi: 10.5281/zenodo.3552717</li>",
            params.genome_only ? "" : "<li>Pertea G, Pertea M. (2020) GFF Utilities: GffRead and GffCompare. F1000Research, 9:304. doi: 10.12688/f1000research.23297.2</li>",
            "<li>Manni M, Berkeley MR, Seppey M, Simão FA, Zdobnov EM. (2021) BUSCO Update: Novel and Streamlined Workflows along with Broader and Deeper Phylogenetic Coverage for Scoring of Eukaryotic, Prokaryotic, and Viral Genomes. Mol Biol Evol, 38(10), 4647-4654. doi: 10.1093/molbev/msab199</li>",
            "<li>Gurevich A, Saveliev V, Vyahhi N, Tesler G. (2013) QUAST: quality assessment tool for genome assemblies. Bioinformatics, 29(8), 1072-1075. doi: 10.1093/bioinformatics/btt086</li>",
            "<li>Rhie A, Walenz BP, Koren S, Phillippy AM. (2020) Merqury: reference-free quality, completeness, and phasing assessment for genome assemblies. Genome Biol, 21(1), 245. doi: 10.1186/s13059-020-02134-9</li>",
            !params.skip_tidk ? "<li>Brown MR, González de la Rosa PM, Blaxter M. (2025) tidk: a toolkit to rapidly identify telomeric repeats from genomic datasets. Bioinformatics, 41(2), btaf049. doi: 10.1093/bioinformatics/btaf049</li>" : "",
            (params.gxdb || params.gxdb_manifest) ? "<li>Astashyn A, Tvedte ES, Sweeney D, et al. (2024) Rapid and sensitive detection of genome contamination at scale with FCS-GX. Genome Biol, 25(1), 60. doi: 10.1186/s13059-024-03198-7</li>" : "",
            (params.gxdb || params.gxdb_manifest) ? "<li>Karlicki M, Antonowicz S, Karnkowska A. (2022) Tiara: deep learning-based classification system for eukaryotic sequences. Bioinformatics, 38(2), 344-350. doi: 10.1093/bioinformatics/btab672</li>" : "",
            params.ortho_version == 'v2' ? "<li>Emms DM, Kelly S. (2019) OrthoFinder: phylogenetic orthology inference for comparative genomics. Genome Biol, 20(1), 238. doi: 10.1186/s13059-019-1832-y</li>" : "<li>Emms DM, Liu Y, Belcher L, et al. (2026) OrthoFinder: improved phylogenetic orthology inference with enhanced accuracy and scalability. Nat Methods. doi: 10.1038/s41592-026-03126-6</li>",
            "<li>Hao Z, Lv D, Ge Y, et al. (2020) RIdeogram: drawing SVG graphics to visualize and map genome-wide data on the idiograms. PeerJ Comput Sci, 6, e251. doi: 10.7717/peerj-cs.251</li>",
            "<li>Yu G, Smith DK, Zhu H, Guan Y, Lam TT. (2017) ggtree: an R package for visualization and annotation of phylogenetic trees with their covariates and other associated data. Methods Ecol Evol, 8(1), 28-36. doi: 10.1111/2041-210X.12628</li>",
            "<li>Paradis E, Schliep K. (2019) ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. Bioinformatics, 35(3), 526-528. doi: 10.1093/bioinformatics/bty633</li>",
            "<li>Wickham H. (2016) ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York.</li>",
            "<li>Lawrence M, Huber W, Pagès H, et al. (2013) Software for Computing and Annotating Genomic Ranges. PLoS Comput Biol, 9(8), e1003118. doi: 10.1371/journal.pcbi.1003118</li>",
            "<li>Wickham H, Averick M, Bryan J, et al. (2019) Welcome to the tidyverse. J Open Source Softw, 4(43), 1686. doi: 10.21105/joss.01686</li>",
            "<li>The pandas development team. pandas-dev/pandas: Pandas. Zenodo. doi: 10.5281/zenodo.3509134</li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
