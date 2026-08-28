# nf-core/genomeqc: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0 - [28 August 2026]

First release of nf-core/genomeqc, which compares the quality of multiple genomes and their annotations. The pipeline runs in two modes depending on the inputs provided: **Genome only** (FASTA files) and **Genome and Annotation** (FASTA plus GTF/GFF files).

### `Added`

- Input genomes and annotations from local files or NCBI accessions, downloaded with [NCBI genome download](https://github.com/kblin/ncbi-genome-download).
- Genome completeness assessment with [BUSCO](https://busco.ezlab.org/), including a **BUSCO Ideogram** plotting the location of markers on the assembly.
- Read-based genome completeness evaluation with [Merqury](https://github.com/marbl/merqury) (optional).
- Telomeric repeat identification and visualisation with [tidk](https://github.com/tolkit/telomeric-identifier) (optional).
- Assembly contiguity and integrity statistics (N50, N90, GC%, number of sequences) with [QUAST](https://github.com/ablab/quast).
- Contamination screening with [FCS-GX](https://github.com/ncbi/fcs/wiki/FCS-GX-quickstart), [FCS-adaptor](https://github.com/ncbi/fcs/wiki/FCS-adaptor-quickstart), and [Tiara](https://ibe-uw.github.io/tiara/), shown as contamination stats in the tree summary plot.
- Annotation statistics with [AGAT](https://agat.readthedocs.io/en/latest/), including a **Gene Overlaps** analysis counting overlapping genes.
- Extraction of the longest protein isoform with [GffRead](https://github.com/gpertea/gffread).
- Orthologous gene inference with [OrthoFinder](https://github.com/davidemms/OrthoFinder) (v2 or v3).
- Repeat/transposable element annotation with RepeatModeler and RepeatMasker (or HITE).
- BUSCO marker-based and orthology-based phylogenetic **Tree Summary** plots combining assembly and annotation summary statistics, with multiple layout styles including a circular layout with configurable quality-threshold rings.
- An executable Shiny app for interactively adjusting the tree plot and summary statistics, with PNG/SVG export.
- HTML and Excel summary reports.
- Aggregated quality-control report with [MultiQC](http://multiqc.info).
- Genome-mode and protein-mode BUSCO/validated-GFF outputs are now published to separate paths so a run using both modes no longer has one mode's files silently overwrite the other's.

### `Fixed`

- Fixed BUSCO not showing in the tree plot.
- Fixed the Shiny app launcher pulling its container from Docker Hub instead of quay.io, which caused a `pull access denied` error.
- Circular layout: quality statistics (sequence count, N50, BUSCO complete/duplicated) are now scored as a colour-vision-safe traffic light (Good/Warn/Poor) against phylogenetic-group thresholds (`--quality_preset`), rather than a sequential colour ramp that misleadingly implied "dark = good". `--show_ring_values` can print each value on its ring.
- Renamed `ORTHOLOGOUS_CHROMOSOMES` to `ORTHO_SEQ_COUNT` ([#196](https://github.com/nf-core/genomeqc/issues/196)), and its outputs accordingly, since the input assembly need not be chromosome-level.
- Fixed `NCBIGENOMEDOWNLOAD` failing under singularity/apptainer with a container-image error.
- Fixed `ORTHOFINDER` silently succeeding when it failed to produce `Orthogroups/Orthogroups.tsv`.
- Fixed `ORTHO_SEQ_COUNT` mapping zero genes for GFFs using a `transcript` feature instead of `mRNA` (common in AUGUSTUS output) ([#174](https://github.com/nf-core/genomeqc/issues/174)).
- Fixed `ORTHOFINDERV2` crashing under `-profile conda` with a misleading "out of RAM" message.
- Fixed `-profile conda` failing to resolve for `RM_DOWNLOAD_DB` and `REPEATMODELER_REPEATMODELER` (HITE remains docker/singularity-only under conda).
- Fixed `HTML_REPORT`/`EXCEL_REPORT` crashing on samplesheets that mix genome-only and annotated assemblies.
- Fixed `gene_overlaps.R` silently dropping genes with unresolved strand, and never detecting overlaps between genes on opposite strands, undercounting overlap statistics.
- Fixed AGAT `sp_statistics` results being effectively absent from both reports; added a feature-first "Annotation stats" tab/sheet and an "Assembly stats" tab to the HTML report.
- Fixed a stale filename reference in the Shiny app that silently broke the ortho-seqs panel.
- Fixed the circular tree layout's rings (`--tree_style circular`) rendering with wildly uneven thickness, and `--show_ring_values`' printed values landing on the wrong ring.
- Fixed the conventional tree layout's TE column title rendering below its panel instead of above it.
- Fixed the FCS-GX panel/ring being drawn against the wrong species, and occasionally showing a fabricated contamination grade for species with no FCS-GX data, whenever FCS-GX only ran for some samples.
- Fixed SVG export failing in the interactive Shiny app.
- Fixed `--RM_download_db true` never actually downloading a DFAM partition, so the parameter had no effect.
- Fixed `GENOMEANNOTATIONBUSCOIDEOGRAM` overwriting outputs when a species is assessed against more than one BUSCO lineage (`--busco_lineage auto`).
- Fixed Merqury outputs (`.qv`, `.completeness.stats`, spectra plots) never being published to `results/merqury/<species>/`.

### `Dependencies`

### `Deprecated`
