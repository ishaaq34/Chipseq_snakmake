# 🐖 Sus scrofa (Pig) ChIP-seq Analysis Pipeline

> **A robust, automated, and interactive Snakemake workflow for ChIP-seq analysis.**  
> *Developed for scalable processing, from raw reads to publication-quality plots.*

---

## 🚀 Quick Start in 3 Steps

If your environment is already set up (Conda/Snakemake installed):

```bash
# 1. Activate your snakemake environment
conda activate snakemake

# 2. Check the configuration (Dry Run)
./dry_run.sh

# 3. Run the full pipeline
./run_pipeline.sh
```

---

## 📖 Overview

This pipeline automates the analysis of ChIP-seq data, specifically tailored for the *Sus scrofa* genome but adaptable to any organism. It handles the entire lifecycle of the data:

1. **Preprocessing**: Adapter trimming and quality filtering.
2. **Alignment**: Mapping reads to the reference genome (Bowtie2).
3. **Refinement**: Deduplication and filtering of blacklisted regions.
4. **Quality Control**: Extensive QC metrics (DeepTools PCA, Correlation, Fingerprints, SPP/IDR).
5. **Peak Calling**: Calling peaks using MACS3 (Narrow & Broad).
6. **Quantification**: Generating BigWigs and signal profiles.
7. **Interactive Exploration**: A custom D3.js dashboard to explore the workflow and results.

---

## 🛠 Prerequisites

You need a Linux/Unix system (or macOS) with **Conda** (or Mamba) installed.

### 1. Install Snakemake

It is recommended to use `mamba` for faster installation:

```bash
conda install -n base -c conda-forge mamba
mamba create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake
```

### 2. External Dependencies

The pipeline builds its own environments for most tools, but you need:

* **Python 3.8+**
* **Graphviz** (for workflow visualization):

    ```bash
    conda install -c conda-forge python-graphviz
    ```

---

## 📂 Input Preparation (The "Hard" Part)

Before running, you must set up your input files correctly.

### 1. Metadata Sheet (`metadata.tsv`)

The pipeline runs entirely based on this file. It defines your experimental design.
**Location**: Defined in `config.yaml` (default: `metadata_from_new_names.tsv`).

**Columns:**

| Column Name | Description | Example |
| :--- | :--- | :--- |
| `raw_id` | **Base name** of your FASTQ file (without `.fastq.gz`). This MUST match the filenames in `raw_reads/`. | `CTRL_R1_IP_H3K9AC` |
| `condition` | Experimental condition (Control, Treated, Infected). | `Control` |
| `replicate` | Biological replicate number. | `1` |
| `assay` | The protein/antibody targeted (H3K9ac, H3K27me3, Input). | `H3K9ac` |
| `mark_type` | Use **"both"** or **"narrow"** / **"broad"** for ChIP samples. Leave **EMPTY** for Input/Control samples. | `both` (or leave empty) |

### 2. Configuration (`config_PE_Sus_scrofa.yaml`)

This YAML file controls all pipeline settings.

* **Genome**: Paths to FASTA/GTF files (`genome:` section).
* **Directories**: Where inputs live and outputs go (`directories:` section).
* **Parameters**: Tool-specific settings for Bowtie2, MACS3, etc.
* **Resources**: Memory allocation for heavy rules.
* **Resources**: Memory allocation for heavy rules.

> **💡 Pro Tip: Automatic Genome Download**
> You don't need to manually download the genome! If the specified `fasta` or `gtf` files are missing, the pipeline will automatically download and decompress them for you.
>
> 1. Set the desired filenames in `fasta: ...` and `gtf: ...`.
> 2. Provide the download links in `fasta_url: ...` and `gtf_url: ...`.

> **⚡️ Efficiency Note: Genome Indexing**
> The pipeline automatically checks for existing Bowtie2 index files (`index/organism*.bt2`). If they are already present, the time-consuming indexing step is **skipped** entirely.

### 3. Raw Data

Place your raw FASTQ files in the directory specified by `raw_reads` in the config (default: `raw_reads/`).

* **Naming Convention**:
  * **Paired-End**: `{raw_id}_R1.fastq.gz` and `{raw_id}_R2.fastq.gz`
  * **Single-End**: `{raw_id}.fastq.gz`

---

## 🗺️ Pipeline Architecture

The workflow executes the following logical steps:

1. **`trim_galore_pe`**: Auto-detection and removal of adapters + quality trimming.
2. **`bowtie2_align`**: Global alignment to *Sus scrofa* genome (very-sensitive mode).
3. **`samtools_sort` & `markdup`**: sorting and removing PCR duplicates (critical for ChIP).
4. **`bam_stats`**: Generates mapping stats, flagstats, and PBC (PCR Bottleneck Coefficient).
5. **`bamCoverage`**: Creates normalized BigWig tracks (RPGC) for visualization in IGV.
6. **`macs3_callpeak`**: Calls peaks relative to Input controls (supporting both narrow and broad peaks).
7. **`DeepTools QC`**:
    * **PCA**: Principal Component Analysis of sample variance.
    * **Correlation**: Heatmap of sample similarity.
    * **Fingerprint**: Quality check for enrichment signal.
    * **Profiles**: Metagene plots (TSS and Gene Body).
8. **`IDR` / `SPP`**: (Optional) Reproducibility analysis between replicates.

---

## 📊 Outputs & Results

The pipeline organizes results into numbered directories for easy navigation:

| Directory | Content |
| :--- | :--- |
| `00_trimmed_reads/` | Cleaned FASTQ files and trimming reports. |
| `01_aligned_bam/` | Sorted BAM files (intermediate). |
| `02_dedup_bam/` | **Final ready-to-use BAM files** (duplicates removed). |
| `03_qc_metrics/` | Text reports (flagstat, PCR metrics) & MultiQC. |
| `04_deeptools_output/` | **BigWig tracks** (`.bw`) and **QC Plots** (PCA, Heatmaps). |
| `06_peak_calling/` | **MACS3 Results** (`.narrowPeak`, `.broadPeak`, `.xls`). |
| `07_frip_scores/` | Fraction of Reads in Peaks (FRiP) calculations. |
| `analysis_summary.xlsx` | **Master Excel Sheet** containing paths to all generated files and QC stats. |
| `docs/results/latest_run.md` | Auto-generated markdown report linking to all plots. |

---

## 🖥️ Interactive Dashboard (New!)

We have included a custom interactive explorer to view the workflow graph and browse results.

**To launch the explorer:**

```bash
python3 generate_interactive_graphs.py --serve
```

(This script will automatically generate the required graph files if they are missing.)

Then open your browser to: **[http://localhost:8000/workflow_explorer.html](http://localhost:8000/workflow_explorer.html)**

* **Workflow Tab**: Zoom and pan through the pipeline architecture.
* **Results Tab**: Browse and view generated plots (PDFs and Interactive HTMLs) directly in the browser.

---

## ⚠️ Troubleshooting

* **"Missing input file"**: Check that your `metadata.tsv` `raw_id` EXACTLY matches your FASTQ filenames in `raw_reads/`. Use `ls raw_reads/` to verify.
* **"Memory Error"**: Increase the `mem_mb` values in `config_PE_Sus_scrofa.yaml` (e.g., `align_mem_mb: 16000`).
* **"Locked Directory"**: If a run crashes, run `snakemake --unlock --configfile config_PE_Sus_scrofa.yaml`.

---
*Generated by Antigravity Agent for automated, reproducible science.*
