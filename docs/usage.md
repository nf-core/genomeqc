# nf-core/genomeqc: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/genomeqc/usage](https://nf-co.re/genomeqc/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

**nf-core/genomeqc** is a pipeline build to aid in the diagnosis of the quality of genome assemblies. It can input several genomes and/or their annotations, and generates metrics such as completeness, contiguity, GC% or number of overlapping genes, which can be later used to assess their quality. Additionally, it outputs a phylogenetic tree plot with the summary metrics. The tree building method uses orthologous genes/BUSCO markers for a quick comparision of metrics across species/samples. This pipeline should not be used for phylogenetic inference.

## Samplesheet input

Before running the pipeline, you will need to create a samplesheet with information about the assemblies you want to process. It has to be a comma-separated file with up to 6 columns, and a header row as shown in the examples below. Use this parameter to specify its location:

```bash
--input '[path to samplesheet file]'
```

The pipeline can be run using NCBI accessions (RefSeq or GenBank) or local files. It needs at least a **genome** (GenBank accession or local file in FASTA format) per assembly to run. If annotations (RefSeq accession or local file in GTF/GFF format) are added, the pipeline will run on both **genomes and annotations**. If you provide local FASTQ reads for an assembly, the pipeline will run Merqury to evaluate completness (see [Running Merqury](#running-merqury)). If you provide taxids in NCBI format and an FCS-GX database, it will run the decontamination subworkflow (see [Running the Decontamination subworkflow](#running-the-decontamination-subworkflow)).

If running the pipeline using **local** files, point to the location these files using the `fasta` and/or `gff` fields:

```csv title="samplesheet.csv"
assembly,fasta,gff
assembly_1,/path/to/genome1.fasta,/path/to/annotation1.gxf
assembly_2,/path/to/genome2.fasta,/path/to/annotation2.gxf
assembly_3,/path/to/genome3.fasta,/path/to/annotation3.gxf
```

If running the pipeline using **ncbi accessions (GenBank and/or RefSeq)**, indicate the corresponding ID using the `ncbi` field:

```csv title="samplesheet.csv"
assembly,ncbi
assembly_1,GCF_000000001.1
assembly_2,GCF_000000002.1
assembly_3,GCF_000000003.1
```

<!-- If running with **Merqury**, you must point to the location of fastq files using the **fastq** field:

```csv title="samplesheet.csv"
assembly,fasta,gxf,fastq
assembly_1,/path/to/genome.fasta,/path/to/annotation.gxf,/path/to/reads.fastq
assembly_2,/path/to/genome.fasta,/path/to/annotation.gxf,/path/to/reads.fastq
assembly_3,/path/to/genome.fasta,/path/to/annotation.gxf,/path/to/reads.fastq
```
-->

This is what the complete samplesheet would look like if using NCBI accessions, and adding FASTQ paths and tax IDs:

```csv title="samplesheet.csv"
assembly,fasta,gff,fastq,taxid
assembly_1,/path/to/genome.fasta,/path/to/annotation.gxf,/path/to/reads.fastq,1234
assembly_2,/path/to/genome.fasta,/path/to/annotation.gxf,/path/to/reads.fastq,1245
assembly_3,/path/to/genome.fasta,/path/to/annotation.gxf,/path/to/reads.fastq,4321
```

This is what the complete samplesheet would look like if using local files, and adding FASTQ paths and tax IDs:

```csv title="samplesheet.csv"
assembly,ncbi,fastq,taxid
assembly_1,GCF_000000001.1,/path/to/reads.fastq,1234
assembly_2,GCF_000000002.1,/path/to/reads.fastq,1245
assembly_3,GCF_000000003.1,/path/to/reads.fastq,4321
```

You can mix different input types in the same samplesheet. If a specific field doesn’t apply to a row, leave it empty (as shown below). The pipeline will automatically detect the input type for each row and run the appropiate subworkflow:

```csv title="samplesheet.csv"
assembly,ncbi,fasta,gff,fastq,taxid
assembly_1,,/path/to/genome.fasta,/path/to/annotation.gxf,/path/to/reads.fastq,1234
assembly_2,,/path/to/genome.fasta,/path/to/annotation.gxf,,4321
assembly_3,,/path/to/genome.fasta,/path/to/annotation.gxf,
assembly_4,,/path/to/genome.fasta,,/path/to/reads.fastq,1245
assembly_5,,/path/to/genome.fasta,,
assembly_6,,/path/to/genome.fassta,,
assembly_7,GCF_000000007.1,,,/path/to/reads.fastq
assembly_8,GCF_000000008.1,,,
assembly_9,GCA_000000009.1,,,/path/to/reads.fastq
assembly_10,GCA_000000010.1,,,,1324
```

As for now, the pipeline doesn't support SRA accession for **Merqury**. We will consider this option the future.

| Column     | Description                                                                                                                         |
| ---------- | ----------------------------------------------------------------------------------------------------------------------------------- |
| `assembly` | Assembly name or custom sample name. Spaces in sample names are automatically converted to underscores (`_`).                       |
| `ncbi`     | ncbi acession. Can be GenBank (starts with `GCA`) or RefSeq (starts with `GCF`).                                                    |
| `fasta`    | Full path to the genome FASTA file. Can be compressed or uncompressed.                                                              |
| `gff`      | Full path to the genome annotation GFF/GTF file. Can be compressed or uncompressed.                                                 |
| `fastq`    | Full path to FASTQ file for long reads (e.g. PacBio or ONT). File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz". |
| `taxid`    | Species taxid for decontamination screening, must be a valid NCBI taxid (numeric string without spaces).                            |

An [example samplesheet](./assets/samplesheet.csv) is provided with the pipeline.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/genomeqc --input ./samplesheet.csv --outdir ./results -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/running/run-pipelines#configuring-pipelines), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/genomeqc -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Modes

#### Genome only

This is the minimal run. The pipeline will run on genome only mode if these inputs are provided in the samplesheet:

1. Path to `fasta` OR
2. `ncbi` GenBank accession.

The pipeline will produce a tree plot summary that can be modified in real time using a packaged shiny app, as well as a MultiQC report with quality statistics.

#### Genome and annotation

The pipeline will run on genome and annotation mode if these inputs are provided in the samplesheet:

1. Path to `fasta` AND path to `gff` OR
2. `ncbi` RefSeq accession.

The pipeline will produce a tree plot summary that can be modified in real time using a packaged shiny app, as well as a MultiQC report with quality statistics.

### Running tidk

The pipeline will run tidk by default. Users can provide a DNA string correspoding to the telomeric repeat of the assemblies of interest using the `--repeat` flag. If no DNA string is provided, tidk will explore the genome and try to compute the telomeric repeat in its cannonical form. E.g.:

```bash
nextflow run nf-core/genomeqc \
  --input ./samplesheet.csv \
  --outdir ./results \
  --repeat ATCG
  -profile docker
```

Users can skip tidk using the `--skip_tidk` flag.

Refer to the [GitHub page](https://github.com/tolkit/telomeric-identifier) for more information about tidk.

### Running Merqury

Users can also run the pipeline using Merqury by supplying the path to sequencing reads under the `fastq` field. Merqury needs both `fasta` and the corresponding `fastq` to run (one single `fastq` per assembly). For paired-end reads, merge/concatenate R1 and R2 into a single file before providing it.

The k-mer size used to build the Meryl database can be tuned with `--kvalue` (default: `21`). For highly heterozygous or large genomes, a larger k (e.g. `31`) may give better completeness estimates.

Refer the [GitHub page](https://github.com/marbl/merqury) for more information on Merqury.

### Running Decontamination

If an NCBI taxid is provided for an assembly in the samplesheet through the `taxid` field, and the path to the FCS-GX database or a manifest to download and build it is given via `--gxdb` or `--gxdb_manifest` respectively, the pipeline will run the decontamination subworkflow. E.g.:

```bash
nextflow run nf-core/genomeqc \
  --input ./samplesheet.csv \
  --outdir ./results \
  --gxdb /patho/to/fcs-gx/all
  -profile docker
```

or

```bash
nextflow run nf-core/genomeqc \
  --input ./samplesheet.csv \
  --outdir ./results \
  --gxdb_manifest 'https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/latest/latest.manifest'
  -profile docker
```

> [!WARNING]
> **The FCS-GX database is considerable in size (~500 GB).** Downloading it in advance over using the manifest is highly recommended.

The decontamination subsworkflow consists of three modules:

- [FCS-GX](https://github.com/ncbi/fcs/wiki/FCS-GX-quickstart): Detection and removal of foreign organisms contamination. Requires the FCS-GX database.
- [FCS-adaptor](https://github.com/ncbi/fcs/wiki/FCS-adaptor-quickstart): Detection and removal of adaptor and vector contamination.
- [Tiara](https://ibe-uw.github.io/tiara/): For DNA sequence classification in two stages:
  1. The sequences are classified to either archaea, bacteria, prokarya, eukarya, organelle or unknown.
  2. The sequences labeled as organelle in the first stage are classified to either mitochondria, plastid or unknown.

> [!WARNING]
> **FCS-GX requires a large database (~500 GB) and is extremely slow when reading it from a standard disk.**
>
> To improve preformance, use `--ramdisk` to point to a `tmpfs` or `ramfs` mount (e.g. `/dev/shm` on Linux). The pipeline will copy the database into memory before each run and clean it up afterwards

When decontamination is run, each genome's FCS-GX contamination percentage is also added to the [tree summary plot](#tree-summary) as a pie chart (conventional layout) or a Good/Warn/Poor traffic-light ring (`--tree_style circular`) - see [Circular tree quality scoring](#circular-tree-quality-scoring).

### Running TE annotation

The pipeline supports two optional methods for transposable element (TE) annotation, selected with the `--te` parameter. TE annotation is skipped by default. Enabling either method adds a TE composition bar to the [tree summary plot](output.md#tree-summary) for each genome.

#### `--te hite`

Runs [HiTE](https://github.com/BioinformaticsToolsmith/HiTE), a fast alignment-free TE identification and masking tool.

> HiTE does not support `-profile conda`/`-profile mamba` - use `-profile docker`, `-profile singularity`, or `-profile podman` instead.

```bash
nextflow run nf-core/genomeqc \
   --input samplesheet.csv \
   --outdir results \
   --te hite \
   -profile docker
```

For plant genomes, also pass `--is_plant true`:

```bash
nextflow run nf-core/genomeqc \
   --input samplesheet.csv \
   --outdir results \
   --te hite \
   --is_plant true \
   -profile docker
```

#### `--te repeatmasker`

Runs a curated TE masking pipeline using the [DFAM](https://www.dfam.org) repeat library:

1. **famdb.py** – extracts a curated repeat library from DFAM h5 partition files (downloaded automatically by default).
2. **Clustering** – deduplicates the library using MMseqs2 `easy-linclust` by default (see [Clustering options](#repeat-library-clustering-options) below).
3. **RepeatMasker** – masks the genome. Runs in rush mode (`-qq`) by default for speed (see [RepeatMasker speed](#repeatmasker-speed) below).

```bash
nextflow run nf-core/genomeqc \
   --input samplesheet.csv \
   --outdir results \
   --te repeatmasker \
   -profile docker
```

To restrict the curated library to a specific taxonomic lineage (strongly recommended — speeds up both the extraction and the masking step), use `--famdb_lineage`:

```bash
nextflow run nf-core/genomeqc \
   --input samplesheet.csv \
   --outdir results \
   --te repeatmasker \
   --famdb_lineage hymenoptera \
   -profile docker
```

By default, the pipeline downloads partition `0` of the DFAM full database. To download additional partitions (required for full taxonomic coverage beyond the root lineage), pass them via `--RM_db`:

```bash
--RM_db "['https://www.dfam.org/releases/current/families/FamDB/dfam39_full.0.h5.gz',\
           'https://www.dfam.org/releases/current/families/FamDB/dfam39_full.1.h5.gz']"
```

If you already have DFAM h5 partition files on disk (highly recommended), use `--famdb_library` to point at them. The download step is skipped automatically.

`--famdb_library` expects **decompressed `.h5` files** — `.h5.gz` files must be extracted first. Only decompress the partitions relevant to your lineage; decompressing the entire database is unnecessary and expensive (each partition is 10–30 GB uncompressed).

Partition 0 always contains the taxonomy metadata and is required. For most lineages, partitions 0 and 1 are sufficient — check the [DFAM partition guide](https://www.dfam.org/releases/current/families/FamDB/) if you need to identify which partition covers your clade.

```bash
# Decompress only the partitions you need (once, on the cluster)
gunzip -k /path/to/dfam38-1_full.0.h5.gz
gunzip -k /path/to/dfam38-1_full.1.h5.gz
```

Then pass them via a glob:

```bash
nextflow run nf-core/genomeqc \
   --input samplesheet.csv \
   --outdir results \
   --te repeatmasker \
   --famdb_library "/path/to/dfam38-1_full.*.h5" \
   --famdb_lineage hymenoptera \
   -profile docker
```

#### Adding de novo repeat discovery with RepeatModeler

By default, `--te repeatmasker` only uses the curated DFAM library. To also run [RepeatModeler](https://www.repeatmasker.org/RepeatModeler/) for de novo discovery, add `--run_repeatmodeler`:

```bash
nextflow run nf-core/genomeqc \
   --input samplesheet.csv \
   --outdir results \
   --te repeatmasker \
   --run_repeatmodeler \
   --famdb_lineage hymenoptera \
   -profile docker
```

The RepeatModeler de novo library is merged with the famdb curated library before masking, giving broader coverage at the cost of significant runtime.

> [!WARNING]
> RepeatModeler is slow — it typically requires 24 CPUs and 24–48 hours per genome. Only enable it if you need de novo discovery beyond the curated DFAM families for your lineage.

#### Repeat library clustering options

Before RepeatMasker runs, the repeat library is deduplicated by a clustering step. Three tools are available via `--te_clusterer`:

| Value      | Tool                    | Notes                                 |
| ---------- | ----------------------- | ------------------------------------- |
| `linclust` | MMseqs2 `easy-linclust` | **Default.** Linear-time, fastest.    |
| `mmseqs`   | MMseqs2 `easy-cluster`  | Slower than linclust, more sensitive. |
| `cdhit`    | CD-HIT-EST              | Traditional approach.                 |

Two thresholds can be tuned:

- `--te_cluster_identity` – minimum sequence identity (default `0.8`). Passed as `-c` to CD-HIT-EST and `--min-seq-id` to MMseqs2.
- `--te_cluster_coverage` – minimum alignment coverage of the shorter sequence (default `0.8`). Passed as `-aS` to CD-HIT-EST and `-c --cov-mode 1` to MMseqs2.

> [!NOTE]
> When running without `--run_repeatmodeler`, clustering runs **once** for the whole lineage library and the result is shared across all genomes. When RepeatModeler is enabled, clustering runs per genome (each genome has a unique de novo library merged in).

#### RepeatMasker speed

RepeatMasker sensitivity can be controlled with `--repeatmasker_speed`:

| Value     | Flag     | Notes                                                               |
| --------- | -------- | ------------------------------------------------------------------- |
| `qq`      | `-qq`    | **Default.** Rush mode — fastest, lowest sensitivity.               |
| `q`       | `-q`     | Quick mode — ~5× faster than default, slightly reduced sensitivity. |
| `default` | _(none)_ | Full sensitivity — slowest.                                         |

For TE quantification in a comparative genomics context, `qq` is usually sufficient. Use `default` if you need a publication-quality masked assembly.

### Circular tree quality scoring

With `--tree_style circular`, sequence count, N50, BUSCO completeness/duplication, and (if decontamination was run) FCS-GX contamination are scored Good/Warn/Poor against `--quality_preset` (default `generic`; also `vertebrate`, `insect`, `plant`, `fungi`, `bacteria`). Pick the preset matching your taxa: e.g. `bacteria` expects a scaffold N50 ≥ 500 kb and ≤ 100 sequences for "Good", while `vertebrate` expects ≥ 10 Mb and ≤ 1000 sequences. These cut-offs are starting points and should be tuned per project. The FCS-GX contamination cut-offs (≥99.5% non-contaminant = Good, 98-99.5% = Warn, <98% = Poor) are fixed across every preset, since there's no taxon-dependent expectation for how much of an assembly should be foreign contamination.

To override individual cut-offs without picking a whole new preset, use `--quality_thresholds` with comma-separated `metric=good:warn` pairs, e.g.:

```bash
--quality_preset bacteria --quality_thresholds 'n50=2e6:5e5,seq_number=50:500'
```

Metrics: `busco_complete`, `busco_duplicated`, `n50`, `seq_number`, `fcs_noncontam`. Direction (higher/lower is better) is fixed per metric, so only the cut-off values are given; any metric not mentioned keeps the selected preset's cut-offs. The same overrides are available interactively via the "Custom..." option in the Shiny app's quality-preset selector.

### The Shiny App

The pipeline outputs an executable that will open a shiny app in your web browser once executed. The app allows the user to change the tree plot parameters in real time, as well as append and remove summary plot statistics. The modified plot can be saved as a png/svg file.

To open the shiny app, do:

```
bash ./results_folder/shiny/app/shiny_app.sh
```

You will need docker to run the shiny app, as it comes packaged in a docker container.

### Running tests

The pipeline can be ran using different test profiles:

1. `-profile test` — Runs genome and annotation with Merqury using **RefSeq accessions** and local FASTQ files.
2. `-profile test_gca` — Runs genome and annotation using **GenBank accessions** (`GCA_`), to confirm the pipeline handles that accession namespace as well as RefSeq (`GCF_`).
3. `-profile test_local` — Runs genome and annotation on local files (`fasta` and `gff`).
4. `-profile test_genomeonly` — Runs genome only on local files (`fasta`).
5. `-profile test_nofastq` — Runs genome and annotation using **RefSeq accessions** (no `fastq`/Merqury).
6. `-profile test_decon` — Tests the decontamination subworkflow (FCS-GX, FCS-Adaptor, Tiara) using the NCBI FCS test-only database.
7. `-profile test_te` — Tests TE annotation with `--te repeatmasker` and `--run_repeatmodeler` on a minimal genome.
8. `-profile test_full` — Runs the full pipeline on a set of Hymenoptera genomes.

Test files are stored in the genomeqc branch of the [test-dataset repository](https://github.com/nf-core/test-datasets/tree/genomeqc).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/genomeqc
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/genomeqc releases page](https://github.com/nf-core/genomeqc/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://charliecloud.io/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#set-max-resources) and [customise process resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#customize-process-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#update-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#modifying-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
