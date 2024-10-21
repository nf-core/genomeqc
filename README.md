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

The pipeline takes a list of genomes and/or annotations (from raw files or Refseq IDs), and runs commonly used tools to assess their quality.

There are different ways you can run this pipeline. 1. Genome only, 2. Annotation only, or 3. Genome and Annotation.

<!-- TODO nf-core:
For an example, see https://github.com/nf-core/rnaseq/blob/master/README.md#introduction
-->

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/contributing/design_guidelines#examples for examples.   
-->

<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline 
Original automatic steps from nf-core pipeline create.
1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Present QC for raw reads ([`MultiQC`](http://multiqc.info/))
-->

**Genome Only:**
1. Downloads the genome files from NCBI `[NCBIGENOMEDOWNLOAD]` - Or you provide your own genomes
2. Describes genome assembly:
2a. `[BUSCO_BUSCO]`: Determines how complete is the genome compared to expected.
2b. `[QUAST]`: Determines the N50, how contiguous the genome is.
2c. More options
3. Summary with MulitQC.

**Genome and Annnotation:**
1. Downloads the genome and gene annotation files from NCBI `[NCBIGENOMEDOWNLOAD]` - Or you provide your own genomes/annotations
2. Describes genome assembly:
2a. `[BUSCO_BUSCO]`: Determines how complete is the genome compared to expected.
2b. `[QUAST]`: Determines the N50, how contiguous the genome is.
2c. More options
3. Describes your annotation : `[AGAT]`: Gene, feature, length, averages, counts. 
4. Extract longest protein fasta sequences `[GFFREAD]`.
5. Finds orthologous genes `[ORTHOFINDER_CAFE]`.
6. Summary with MulitQC.

**Annnotation Only:**
1. Downloads the gene annotation files from NCBI `[NCBIGENOMEDOWNLOAD]` - Or you provide your own annotations.
2. Describes your annotation : `[AGAT]`: Gene, feature, length, averages, counts.
3. Summary with MulitQC.

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a `samplesheet.csv`, where your input data points to genomes + or annotations:

```csv
sample,genome.fasta,annotation.gff
Homo_sapiens,/path/to/genome.fasta,/path/to/annotation.gff3
Gorilla_gorilla,/path/to/genome.fasta,
Pan_paniscus,,/path/to/annotation.gff3
```

Or to Refseq IDs of your species:

```csv
sample,refseqID
Pongo_abelii,GCF_028885655.2
Macaca_mulatta,GCF_003339765.1
```

You can mix the two input types. Also, notice you can leave the genome or annotation absent.

Each row represents a species, with its associated genome, gff or Refseq ID (to autodownload the genome + gff).

Now, you can run the pipeline using:

```
nextflow run main.nf -resume -profile docker,test --outdir results
```

or 

```
nextflow run main.nf -resume -profile docker --input assets/samplesheet.csv --outdir results
```

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

```bash
nextflow run ecoflow/genomeqc \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

## Credits

ecoflow/genomeqc was originally written by Chris Wyatt, Fernando Duarte.

We thank the following people for their extensive assistance in the development of this pipeline:

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
