ChIP-seq Analysis Workflow
===========================

This workflow performs comprehensive ChIP-seq analysis from raw FASTQ files to consensus peaks and motif discovery.

Features
--------

* Quality control and trimming with Trim_galore
* Alignment with Bowtie2 and duplicate removal
* BigWig generation with configurable normalization
* Peak calling with MACS3 (narrow and broad peaks)
* Reproducibility analysis with IDR
* Consensus peak calling (high-confidence peaks)
* Motif discovery with HOMER (optional)
* Condition-specific profile plots (TSS and gene body)

Configuration
-------------

The workflow is configured via ``config.yaml`` and uses a metadata TSV file for sample definitions.

Report Contents
---------------

This report includes:

* Quality control metrics
* Peak calling results
* Reproducibility analysis (IDR)
* Consensus peaks per condition
* Visualization plots (correlation, PCA, profiles)
* Motif discovery results (if enabled)
