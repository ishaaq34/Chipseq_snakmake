# ChIP-seq Analysis Workflow

A comprehensive Snakemake workflow for ChIP-seq data analysis from raw FASTQ files to publication-ready results.

## Features

- ✅ **Flexible sequencing support**: Single-end and paired-end
- ✅ **Metadata-driven**: Define experiments in simple TSV format
- ✅ **Automatic genome setup**: Optional download of FASTA/GTF from URLs
- ✅ **Quality control**: Trim_galore, MultiQC, DeepTools QC plots
- ✅ **Peak calling**: MACS3 with narrow/broad peak support
- ✅ **Reproducibility analysis**: IDR and consensus peaks
- ✅ **Visualization**: BigWig generation, TSS/gene body profiles
- ✅ **Interactive reports**: Datavzrd tables and automatic HTML reports
- ✅ **Motif discovery**: HOMER integration (optional)
- ✅ **Fully configurable**: All parameters in `config.yaml`

## Quick Start

### 1. Setup

```bash
# Clone or copy this workflow
cd chipseq_workflow

# Install dependencies (if not using conda)
conda env create -f environment.yaml
conda activate chipseq
```

### 2. Configure Your Analysis

Edit `config.yaml`:

```yaml
# Reference genome
genome:
  fasta: "genome.fa"
  gtf: "genes.gtf"
  index_prefix: "genome_index/genome"
  effective_size: 2700000000
  
  # Optional: auto-download if URLs provided
  fasta_url: ""
  gtf_url: ""

# Sequencing type: "paired" or "single"
sequencing_type: "paired"

# Metadata TSV (defines your experiment)
macs3:
  metadata_tsv: "samples_metadata.tsv"
```

### 3. Create Metadata TSV

`samples_metadata.tsv`:

```tsv
raw_id condition replicate assay mark_type
H3K9ac_rep1 treated 1 H3K9ac narrow
H3K9ac_rep2 treated 2 H3K9ac narrow
Input_rep1 treated 1 input 
Input_rep2 treated 2 input 
```

**Required columns:**

- `raw_id`: Sample identifier (must match FASTQ filenames)
- `condition`: Experimental condition
- `replicate`: Replicate number (1, 2, 3...)
- `assay`: ChIP mark or "input"
- `mark_type`: "narrow" or "broad" (for ChIP samples)

### 4. Place FASTQ Files

**Single-end:**

```
Raw_reads/
├── H3K9ac_rep1.fastq.gz
├── H3K9ac_rep2.fastq.gz
├── Input_rep1.fastq.gz
└── Input_rep2.fastq.gz
```

**Paired-end:**

```
Raw_reads/
├── H3K9ac_rep1_R1.fastq.gz
├── H3K9ac_rep1_R2.fastq.gz
├── H3K9ac_rep2_R1.fastq.gz
├── H3K9ac_rep2_R2.fastq.gz
└── ...
```

### 5. Run Workflow

```bash
# Dry run to check
snakemake -n

# Run with 8 cores
snakemake --cores 8

# Run with conda environments
snakemake --use-conda --cores 8

# Generate final report
snakemake --report chipseq_report.zip
```

## Output Structure

```
chipseq_workflow/
├── Trim_galore/              # Trimmed FASTQ files
├── aligned_bam/              # Sorted BAM files
├── de_duplicate4/            # Deduplicated BAM files
├── bam_QC_PBC/              # QC metrics + MultiQC report
├── deeptools_QC/            # BigWigs, correlation, PCA
├── macs3_results/           # Peak calls per replicate
├── frip_analysis/           # FRiP metrics + interactive table
├── idr_analysis/            # IDR results
├── consensus_peaks/         # High-confidence reproducible peaks
├── bigwig_averaged/         # Mean signal per condition
├── bigwig_normalized/       # log2(ChIP/Input) tracks
├── condition_profiles/      # TSS & gene body plots
└── homer_motifs/            # Motif discovery (optional)
```

## Key Outputs

### Quality Control

- `bam_QC_PBC/alignment_qc_report.html` - MultiQC report
- `deeptools_QC/correlation_heatmap.pdf` - Sample correlation
- `deeptools_QC/pca.pdf` - PCA plot
- `deeptools_QC/fingerprints.pdf` - Signal enrichment

### Peak Calling

- `macs3_results/{sample}_peaks.narrowPeak` - Per-replicate peaks
- `frip_analysis/interactive_table/index.html` - FRiP metrics (interactive)
- `consensus_peaks/{condition}_{assay}_consensus.bed` - Reproducible peaks

### Visualization

- `bigwig_normalized/{condition}_{assay}_log2_over_input.bw` - Normalized tracks
- `condition_profiles/{condition}_{assay}_TSS_log2_profile.pdf` - TSS enrichment

### Final Report

- `chipseq_report.zip` - Self-contained HTML with all results

## Configuration Options

### Tool-Specific Parameters

**Trim_galore:**

```yaml
trim_galore:
  quality: 20
  min_length: 20
  adapter: "auto"
```

**Bowtie2:**

```yaml
bowtie2:
  sensitivity: "--very-sensitive"
  no_unal: true
  no_mixed: true       # PE only
  no_discordant: true  # PE only
```

**Samtools:**

```yaml
samtools:
  mapq_threshold: 30
  sam_flags_exclude: 768  # 256 (secondary) + 512 (QC-fail)
```

**MACS3:**

```yaml
macs3:
  params:
    genome_size: "hs"
    pvalue_cutoff: 0.01
    keep_dup: "all"
    broad_cutoff: 0.1
```

**DeepTools:**

```yaml
deeptools:
  bin_size: 1000
  normalization_method: "RPGC"
  tss_upstream: 3000
  tss_downstream: 10000
```

## Advanced Usage

### Custom Mark Types

Edit metadata `mark_type` column:

- `narrow`: Sharp peaks (H3K4me3, H3K9ac, H3K27ac, CTCF)
- `broad`: Diffuse peaks (H3K27me3, H3K36me3, H3K9me3)

### Enable Optional Tools

```yaml
optional_tools:
  run_spp: true
  spp_script: "/path/to/run_spp.R"
  
  run_homer: true
  homer_params:
    size: 200
    length: "8,10,12"
```

### Generate Reports

```bash
# Full HTML report
snakemake --report chipseq_report.zip

# Unzip and view
unzip chipseq_report.zip
open report.html
```

## Troubleshooting

**"No input found for sample X"**
→ Ensure each ChIP sample has matching input with same condition + replicate

**"Unbalanced replicates"**
→ All assays in a condition must have same number of replicates

**Dry run fails with wildcard error**
→ Normal for `build_index` rule, workflow will run properly

## Citation

If you use this workflow, please cite:

- Snakemake: Mölder et al., F1000Research 2021
- MACS3: Zhang et al., Genome Biology 2008
- DeepTools: Ramírez et al., Nucleic Acids Research 2016
- Datavzrd: Köster et al., 2023

## License

MIT License - See LICENSE file for details
