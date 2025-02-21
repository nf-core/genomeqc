# ![nf-core/genomeqc](docs/images/nf-core-genomeqc_logo_light.png#gh-light-mode-only) ![nf-core/genomeqc](docs/images/nf-core-genomeqc_logo_dark.png#gh-dark-mode-only)

[![GitHub Actions CI Status](https://github.com/ecoflow/genomeqc/actions/workflows/ci.yml/badge.svg)](https://github.com/ecoflow/genomeqc/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/ecoflow/genomeqc/actions/workflows/linting.yml/badge.svg)](https://github.com/ecoflow/genomeqc/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/ecoflow/genomeqc)

## Introduction

**ecoflow/genomeqc** is a bioinformatics pipeline that compares the quality of multiple genomes, along with their annotations.

The pipeline takes a list of genomes and annotations (from raw files or Refseq IDs), and runs commonly used tools to assess their quality.

There are three different ways you can run this pipeline:
 1. Genome only
 2. Annotation only
 3. Genome and Annotation.

**Currently, only Genome plus Annotation is functional**

<!-- TODO nf-core:
For an example, see https://github.com/nf-core/rnaseq/blob/master/README.md#introduction
-->

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/contributing/design_guidelines#examples for examples.
-->

![pipeline_diagram](docs/images/nf-core-genomeqc_metro_map_v2.png)

**1. Genome and Annotation:**
1. Downloads the genome and gene annotation files from NCBI: [NCBI genome download](https://github.com/kblin/ncbi-genome-download) - Or you provide your own genomes/annotations
2. Describes genome assembly:
   1. [BUSCO](https://busco.ezlab.org/): Evaluates genome completeness based on **single copy markers**.
   2. **BUSCO Ideogram**: Plots the location of markers on the assembly.
   3. [Merqury](https://github.com/marbl/merqury) (optional): Evaluates genome completeness based on sequencing reads.
   4. [tidk](https://github.com/tolkit/telomeric-identifier) (optional): Indetifies and visualises telomeric repeats.
   5. [QUAST](https://github.com/ablab/quast): Computes contiguity and integrity statistics: N50, N90, GC%, number of sequences.
   6. More options...
3. Describes annotation : 
   1. [AGAT](https://agat.readthedocs.io/en/latest/): Number of genes, features, length...
   2. **Gene Overlaps**: Finds the number of overlapping genes.
   3. More options...
5. Extracts longest protein isoform: [GffRead](https://github.com/gpertea/gffread).
6. Finds orthologous genes: [Orthofinder](https://github.com/davidemms/OrthoFinder).
7. Plots an orthology-based phylogenetic tree : **Tee Summary**, as well as other relevant stats from the above steps.
8. Summary with [MultiQC](http://multiqc.info).

**2. Genome Only (in development):**
1. Downloads the genome files from NCBI: [NCBI genome download](https://github.com/kblin/ncbi-genome-download) - Or you provide your own genomes
2. Describes genome assembly:
   1. [BUSCO](https://busco.ezlab.org/): Evaluates genome completeness based on **single copy markers**.
   2. **BUSCO Ideogram**: Plots the location of markers on the assembly.
   3. [tidk](https://github.com/tolkit/telomeric-identifier) (optional): Indetfies and visualises telomeric repeats.
   3. [QUAST](https://github.com/ablab/quast): Computes contiguity and integrity statistics: N50, N90, GC%, number of sequences.
3. Summary with [MultiQC](http://multiqc.info).

**3. Annnotation Only (in development):**
1. Downloads the gene annotation files from NCBI: [NCBI genome download](https://github.com/kblin/ncbi-genome-download) - Or you provide your own annotations.
2. Describes your annotation : [AGAT](https://agat.readthedocs.io/en/latest/): Gene, feature, length, averages, counts.
3. Summary with [MultiQC](http://multiqc.info).

> [!WARNING]
> We strongly suggest users to specify the lineage using the `--busco_lineage` parameter, as setting the lineage to `auto` (default value) might cause problems with `BUSCO` during the lineage determination step.

> [!NOTE]
> `BUSCO Ideogram` will only plot those chromosomes -or scaffolds- that contain single copy markers.

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare an input **samplesheet** in **csv format** (e.g. `samplesheet.csv`). You can prepare your sampplesheet using:

###  1. Local files

Simply point out to your local genome assembly and annotation (in FASTA and GFF format, respectively) using the `fasta` and `gff` fields, leaving the rest of the fields empty:

```csv
species,refseq,fasta,gff,fastq
Homo_sapiens,,/path/to/genome.fasta,/path/to/annotation.gff3,
Gorilla_gorilla,,/path/to/genome.fasta,/path/to/annotation.gff3,
Pan_paniscus,,/path/to/genome.fasta,/path/to/annotation.gff3,
```

When running on ``--genome_only`` mode, you can leave the `gff` field empty. Otherwise, this field will be ignored.

### 2. Refseq IDs

Additionally, you can run the pipeline using the Refseq IDs of the assemblies of interest in the `refseq` field, leaving the rest of the fields empty:

```csv
species,refseq,fasta,gff,fastq
Pongo_abelii,GCF_028885655.2,,,
Macaca_mulatta,GCF_003339765.1,,,
```

### Run the pipeline

<!--
You can mix the two input types **(in development)**.
-->

Run the pipeline using:

```bash
nextflow run nf-core/genomeqc \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

You can run the pipeline using a test profile and docker:

```bash
nextflow run nf-core/genomeqc -profile docker --outdir ./results
```

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/genomeqc/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/genomeqc/output).

## Credits

ecoflow/genomeqc was originally written by [Chris Wyatt](https://github.com/chriswyatt1) and [Fernando Duarte](https://github.com/FernandoDuarteF) at the University College London.

We thank the following people for their extensive assistance in the development of this pipeline:

- [Mahesh Binzer-Panchal](https://github.com/mahesh-panchal) ([National Bioinformatics Infrastructure Sweden](https://nbis.se/))
- [Usman Rashid](https://github.com/GallVp) ([The New Zealand Institute for Plant and Food Research](https://www.plantandfood.com/en-nz/))
- [Lauren Huet](https://github.com/LaurenHuet) ([Schmidt Ocean Institute](https://schmidtocean.org/))
- [Stephen Turner](https://github.com/stephenturner/) ([Colossal Biosciences](https://colossal.com/))
- [Felipe Perez Cobos](https://github.com/fperezcobos) ([Institute of Agrifood Research and Technology](https://www.irta.cat/en/))

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use ecoflow/genomeqc for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).


python  app/downloader-utility.py --clade "Chordata" --project_name "DToL" --data_status "Mapped Reads - Done" --experiment_type "Chromium genome"  --download_location "/Users/raheela/Documents" --download_option "assemblies" --species_list "Apamea sordens,Bufo bufo"