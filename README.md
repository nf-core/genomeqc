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

There are three different ways you can run this pipeline. 1. Genome only, 2. Annotation only, or 3. Genome and Annotation. **Only Genome plus Annotation is functional**

<!-- TODO nf-core:
For an example, see https://github.com/nf-core/rnaseq/blob/master/README.md#introduction
-->

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/contributing/design_guidelines#examples for examples.   
-->

**Genome and Annnotation:**
1. Downloads the genome and gene annotation files from NCBI `[NCBIGENOMEDOWNLOAD]` - Or you provide your own genomes/annotations
2. Describes genome assembly:
2a. `[BUSCO_BUSCO]`: Determines how complete is the genome compared to expected (protein mode).
2b. `[BUSCO_IDEOGRAM]: Plots the location of BUSCO markers on the assembly.
2c. `[QUAST]`: Determines the N50, how contiguous the genome is.
2d. More options
3. Describes your annotation : `[AGAT]`: Gene, feature, length, averages, counts. 
4. Extract longest protein fasta sequences `[GFFREAD]`.
5. Finds orthologous genes `[ORTHOFINDER]`.
6. Summary with MulitQC.q

> [!WARNING] We strongly suggest users to specify the lineage using the `--busco_lineage` parameter, as setting the lineage to `auto` (value by default) might cause problems with `[BUSCO]` during the leneage determination step.
> [!NOTE] `BUSCO_IDEOGRAM` will plot only those chromosomes -or scaffolds- that contain single copy markers.

**Genome Only (in development):**
1. Downloads the genome files from NCBI `[NCBIGENOMEDOWNLOAD]` - Or you provide your own genomes
2. Describes genome assembly:
2a. `[BUSCO_BUSCO]`: Determines how complete is the genome compared to expected (genome mode).
2b. `[QUAST]`: Determines the N50, how contiguous the genome is.
2c. More options
3. Summary with MulitQC.

**Annnotation Only (in development):**
1. Downloads the gene annotation files from NCBI `[NCBIGENOMEDOWNLOAD]` - Or you provide your own annotations.
2. Describes your annotation : `[AGAT]`: Gene, feature, length, averages, counts.
3. Summary with MulitQC.

In addition to the three different modes described above, it is also possible to run the pipeline with or without sequencing reads. When supplying sequencing reads, Merqury can also be run. [Merqury](https://github.com/marbl/merqury) is a tool for genome quality assessment that uses k-mer counts from raw sequencing data to evaluate the accuracy and completeness of a genome assembly. Meryl is the companion tool that efficiently counts and stores k-mers from sequencing reads, enabling Merqury to estimate metrics like assembly completeness and base accuracy. These tools provide a k-mer-based approach to assess assembly quality, helping to identify potential errors or gaps.​

To run the pipeline with reads, you must supply a single FASTQ file for each genome in the samplesheet, alongside the `--run_merqury` flag. It is assumed that reads used to create the assembly are from long read technology such as PacBio or ONT, and are therefore single end. If reads are in a .bam file, they must be converted to FASTQ format first. If you have paired end reads, these must be interleaved first.

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a `samplesheet.csv`, where your input data points to genomes + or annotations:

```csv
species,refseq,fasta,gff,fastq
Homo_sapiens,,/path/to/genome.fasta,/path/to/annotation.gff3,[/path/to/reads.fq.gz]
Gorilla_gorilla,,/path/to/genome.fasta,/path/to/annotation.gff3,[/path/to/reads.fq.gz]
Pan_paniscus,,/path/to/genome.fasta,/path/to/annotation.gff3,[/path/to/reads.fq.gz]
```

When running on ``--genome_only`` mode, you can leave the **gff** field empty. Otherwise, this field will be ignored.

Additionally, you can run the pipeline using the Refseq IDs of your species:

```csv
species,refseq,fasta,gff,fastq
Pongo_abelii,GCF_028885655.2,,,[/path/to/reads.fq.gz]
Macaca_mulatta,GCF_003339765.1,,,[/path/to/reads.fq.gz]
```

The **fastq** field is optional. Supply sequencing reads if you intend to run merqury using the `--run_merqury`. Otherwise, this filed will be ignored.

You can mix the two input types **(in development)**.

Each row represents a species, with its associated genome, gff or Refseq ID (to autodownload the genome + gff).

You can run the pipeline using test profiles or example input samplesheets. To run a test set with a samplesheet containing reads:

```
nextflow run main.nf -resume -profile docker,test --outdir results --run_merqury
```

To run this pipeline on an example samplesheet included in the repo assets (_does not include reads_):

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

- [Stephen Turner](https://github.com/stephenturner/) ([Colossal Biosciences](https://colossal.com/))

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
