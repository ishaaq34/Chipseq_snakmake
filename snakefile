# =====================================================
# ChIP-seq Analysis Snakemake Pipeline
# =====================================================
# Full workflow: FASTQ → Alignment → Deduplication → QC → Visualization
# Can also start from pre-aligned BAMs (set start_from: "bam" in config)
#
# Usage:
#   snakemake --cores 8 --configfile config.yaml
#   snakemake --cores 8 --configfile config.yaml --use-conda
#
# Requirements:
#   - config.yaml in working directory
#   - sample_id.txt with one sample ID per line
#   - For "fastq" mode: trimmed FASTQ files in trim/ directory
#   - For "bam" mode: sorted BAMs in bowtie_align/ directory
# =====================================================

import os

# =====================================================
# Configuration and Report Setup
# =====================================================
configfile: "config.yaml"

# Global report configuration
report: "report/workflow.rst"

# Extract configuration
GENOME_FASTA = config["genome"]["fasta"]
GTF_FILE = config["genome"]["gtf"]
GENOME_INDEX = config["genome"]["index_prefix"]
EFFECTIVE_GENOME_SIZE = config["genome"]["effective_size"]

RAW_READS_DIR = config["directories"]["raw_reads"]
TRIM_DIR = config["directories"]["trim"]
OUT_ALIGN = config["directories"]["alignment"]
OUT_DEDUP = config["directories"]["dedup"]
OUT_QC = config["directories"]["qc"]
OUT_DEEPTOOLS = config["directories"]["deeptools"]
LOG_DIR = config["directories"]["logs"]

THREADS = config["threads"]
BIN_SIZE = config["deeptools"]["bin_size"]
SMOOTH_LENGTH = config["deeptools"]["smooth_length"]
MAPQ_THRESHOLD = config["parameters"]["mapq_threshold"]
SAM_FLAGS = config["parameters"]["sam_flags_exclude"]
NORMALIZATION_METHOD = config["deeptools"].get("normalization_method", "RPGC")

START_FROM = config.get("start_from", "fastq")
SEQ_TYPE = config.get("sequencing_type", "paired")  # "paired" or "single"

# SPP configuration
SPP_SCRIPT = config["optional_tools"]["spp_script"]
RUN_SPP = config["optional_tools"]["run_spp"] and os.path.exists(SPP_SCRIPT)

# HOMER configuration
RUN_HOMER = config.get("optional_tools", {}).get("run_homer", False)

# =====================================================
# Metadata TSV Validation and Parsing
# =====================================================
def validate_metadata_tsv(tsv_path):
    """Validate metadata table and derive ChIP-input pairs"""
    import pandas as pd
    
    df = pd.read_csv(tsv_path, sep="\t")
    
    # Check required columns
    required = {"raw_id", "condition", "replicate", "assay", "mark_type"}
    if not required.issubset(df.columns):
        raise ValueError(f"Metadata TSV missing columns: {required - set(df.columns)}")
    
    # Check for nulls in key columns
    key_cols = ["raw_id", "condition", "replicate", "assay"]
    if df[key_cols].isnull().any().any():
        raise ValueError("Metadata contains empty cells in required columns")
    
    # Validate assays
    assays = set(df["assay"].unique())
    if "input" not in assays:
        raise ValueError("Metadata must include input control (assay='input')")
    
    chip_assays = assays - {"input"}
    if len(chip_assays) < 1:
        raise ValueError("Metadata must include at least one ChIP assay")
    
    # Check balanced replicates across conditions (per assay)
    for assay in chip_assays:
        assay_df = df[df["assay"] == assay]
        rep_counts = assay_df.groupby("condition")["replicate"].nunique()
        
        if rep_counts.nunique() != 1:
            raise ValueError(f"{assay}: unbalanced replicates across conditions: {rep_counts.to_dict()}")
        
        if rep_counts.iloc[0] < 2:
            raise ValueError(f"{assay}: need ≥2 replicates per condition, found {rep_counts.iloc[0]}")
    
    print(f"✓ Metadata validation passed: {len(df)} rows, {len(assays)} assays, {df['condition'].nunique()} conditions")
    return df

def derive_chip_input_pairs(metadata_df):
    """Automatically pair ChIP with input based on (condition, replicate)"""
    pairs = {}
    
    # Get all ChIP assays
    chip_assays = set(metadata_df[metadata_df["assay"] != "input"]["assay"].unique())
    
    # For each ChIP assay
    for assay in chip_assays:
        assay_df = metadata_df[metadata_df["assay"] == assay]
        
        for _, chip_row in assay_df.iterrows():
            cond = chip_row["condition"]
            rep = chip_row["replicate"]
            
            # Find matching input for this (condition, replicate)
            input_row = metadata_df[
                (metadata_df["condition"] == cond) & 
                (metadata_df["replicate"] == rep) & 
                (metadata_df["assay"] == "input")
            ]
            
            if len(input_row) == 0:
                raise ValueError(f"No input found for {chip_row['raw_id']} (condition={cond}, replicate={rep})")
            
            input_row = input_row.iloc[0]
            
            sample_name = chip_row["raw_id"]
            pairs[sample_name] = {
                "mark_type": chip_row["mark_type"],
                "control": input_row["raw_id"],
                "condition": cond,
                "replicate": rep,
                "assay": assay
            }
    
    # Add input samples (without peak calling)
    for _, input_row in metadata_df[metadata_df["assay"] == "input"].iterrows():
        input_name = input_row["raw_id"]
        if input_name not in pairs:
            pairs[input_name] = {
                "mark_type": None,
                "control": None,
                "condition": input_row["condition"],
                "replicate": input_row["replicate"],
                "assay": "input"
            }
    
    return pairs

# MACS3 configuration - Auto-detect mode
if config.get("metadata_tsv") and os.path.exists(config["metadata_tsv"]):
    # Mode 2: Metadata-driven
    print(f"Using metadata TSV: {config['metadata_tsv']}")
    metadata_df = validate_metadata_tsv(config["metadata_tsv"])
    MACS3_SAMPLES_CONFIG = derive_chip_input_pairs(metadata_df)
    print(f"✓ Auto-paired {len([s for s in MACS3_SAMPLES_CONFIG if MACS3_SAMPLES_CONFIG[s]['mark_type']])} ChIP samples with inputs")
    
    # Extract sample IDs from metadata (no separate samples_file needed!)
    SAMPLES = list(metadata_df["raw_id"].unique())
    print(f"✓ Extracted {len(SAMPLES)} samples from metadata TSV")
    
elif config.get("macs3_samples"):
    # Mode 1: Explicit (backward compatible)
    MACS3_SAMPLES_CONFIG = config["macs3_samples"]
    if MACS3_SAMPLES_CONFIG:
        print(f"✓ Using explicit sample configuration ({len(MACS3_SAMPLES_CONFIG)} samples)")
    
    # Read samples from file (traditional method)
    if "samples_file" in config and os.path.exists(config["samples_file"]):
        SAMPLES = [line.strip() for line in open(config["samples_file"])]
    else:
        raise ValueError("No metadata_tsv or valid samples_file found in config. Please provide one.")
else:
    # No MACS3 configured - must have samples_file
    MACS3_SAMPLES_CONFIG = {}
    if "samples_file" in config and os.path.exists(config["samples_file"]):
        SAMPLES = [line.strip() for line in open(config["samples_file"])]
    else:
        raise ValueError("No metadata_tsv or samples_file found in config. Please provide sample information.")

# Get samples that need peak calling (exclude inputs)
MACS3_SAMPLES = [s for s in MACS3_SAMPLES_CONFIG.keys() 
                 if MACS3_SAMPLES_CONFIG[s].get("mark_type")]

# Helper functions for SE/PE conditional inputs
def get_trim_input(wildcards):
    """Return trim_galore input files based on sequencing type"""
    if SEQ_TYPE == "paired":
        return {
            "r1": os.path.join(RAW_READS_DIR, f"{wildcards.sample}_R1.fastq.gz"),
            "r2": os.path.join(RAW_READS_DIR, f"{wildcards.sample}_R2.fastq.gz")
        }
    else:
        return {"r1": os.path.join(RAW_READS_DIR, f"{wildcards.sample}.fastq.gz")}

def get_trim_output(wildcards):
    """Return trim_galore output files based on sequencing type"""
    if SEQ_TYPE == "paired":
        return {
            "r1": os.path.join(TRIM_DIR, f"{wildcards.sample}_R1_val_1.fq.gz"),
            "r2": os.path.join(TRIM_DIR, f"{wildcards.sample}_R2_val_2.fq.gz"),
            "report1": os.path.join(TRIM_DIR, f"{wildcards.sample}_R1.fastq.gz_trimming_report.txt"),
            "report2": os.path.join(TRIM_DIR, f"{wildcards.sample}_R2.fastq.gz_trimming_report.txt")
        }
    else:
        return {
            "r1": os.path.join(TRIM_DIR, f"{wildcards.sample}_trimmed.fq.gz"),
            "report1": os.path.join(TRIM_DIR, f"{wildcards.sample}.fastq.gz_trimming_report.txt")
        }

def get_align_input(wildcards):
    """Return alignment input files based on sequencing type"""
    if SEQ_TYPE == "paired":
        return {
            "r1": os.path.join(TRIM_DIR, f"{wildcards.sample}_R1_val_1.fq.gz"),
            "r2": os.path.join(TRIM_DIR, f"{wildcards.sample}_R2_val_2.fq.gz")
        }
    else:
        return {"r1": os.path.join(TRIM_DIR, f"{wildcards.sample}_trimmed.fq.gz")}

# Create output directories
os.makedirs(TRIM_DIR, exist_ok=True)
os.makedirs(OUT_ALIGN, exist_ok=True)
os.makedirs(OUT_DEDUP, exist_ok=True)
os.makedirs(OUT_QC, exist_ok=True)
os.makedirs(OUT_DEEPTOOLS, exist_ok=True)
os.makedirs(LOG_DIR, exist_ok=True)

# =====================================================
# Target Outputs
# =====================================================
ALL_BIGWIGS = expand(os.path.join(OUT_DEEPTOOLS, "{sample}.dedup.bw"), sample=SAMPLES)
ALL_FLAGS = expand(os.path.join(OUT_QC, "{sample}.flagstat.txt"), sample=SAMPLES)
ALL_SPP = expand(os.path.join(OUT_QC, "{sample}_spp.qc.txt"), sample=SAMPLES) if RUN_SPP else []

# MultiQC reports
MULTIQC_ALIGNMENT = os.path.join(OUT_QC, "alignment_qc_report.html")

ALL_DEEPTOOLS_QC = [
    os.path.join(OUT_DEEPTOOLS, "all_samples_bins.npz"),
    os.path.join(OUT_DEEPTOOLS, "heatmap_SpearmanCorr.pdf"),
    os.path.join(OUT_DEEPTOOLS, "correlation_matrix.tab"),
    os.path.join(OUT_DEEPTOOLS, "PCA.pdf"),
    os.path.join(OUT_DEEPTOOLS, "fingerprints.pdf"),
    os.path.join(OUT_DEEPTOOLS, "coverage_histogram.pdf")
]

ALL_PROFILE_OUTPUTS = [
    os.path.join(OUT_DEEPTOOLS, "matrix_all_bw_TSS.gz"),
    os.path.join(OUT_DEEPTOOLS, "all_samples_TSS_profile.pdf"),
    os.path.join(OUT_DEEPTOOLS, "matrix_all_bw_scalar.gz"),
    os.path.join(OUT_DEEPTOOLS, "all_samples_scalar_profile.pdf")
]

# MACS3 peak calling outputs (conditional on configuration)
ALL_MACS3_PEAKS = []
if MACS3_SAMPLES:
    for sample in MACS3_SAMPLES:
        mark_type = MACS3_SAMPLES_CONFIG[sample]["mark_type"]
        if mark_type == "narrow":
            ALL_MACS3_PEAKS.extend([
                os.path.join("macs3_results", f"{sample}_peaks.narrowPeak"),
                os.path.join("macs3_results", f"{sample}_summits.bed")
            ])
        elif mark_type == "broad":
            ])

# FRiP calculation outputs
ALL_FRIP = []
if MACS3_SAMPLES:
    for sample in MACS3_SAMPLES:
        ALL_FRIP.append(os.path.join("frip_analysis", f"{sample}.frip.txt"))
    # Aggregate table and interactive view
    ALL_FRIP.append(os.path.join("frip_analysis", "all_samples_frip.tsv"))
    ALL_FRIP.append(os.path.join("frip_analysis", "interactive_table", "index.html"))

# IDR analysis outputs (conditional on metadata TSV with ≥2 replicates)
ALL_IDR = []
ALL_CONSENSUS_PEAKS = []
if MACS3_SAMPLES_CONFIG and config.get("metadata_tsv"):
    # Group samples by condition
    from itertools import combinations
    
    condition_samples = {}
    for sample, info in MACS3_SAMPLES_CONFIG.items():
        if info.get("mark_type"):  # Only ChIP samples
            cond = info.get("condition", "default")
            if cond not in condition_samples:
                condition_samples[cond] = []
            condition_samples[cond].append(sample)
    
    # Generate pairwise IDR comparisons per condition
    for condition, samples in condition_samples.items():
        if len(samples) >= 2:
            for sample1, sample2 in combinations(sorted(samples), 2):
                rep1 = MACS3_SAMPLES_CONFIG[sample1]["replicate"]
                rep2 = MACS3_SAMPLES_CONFIG[sample2]["replicate"]
                
                ALL_IDR.append(os.path.join("idr_analysis", f"{condition}_rep{rep1}_vs_rep{rep2}_idr.txt"))
                ALL_IDR.append(os.path.join("idr_analysis", f"{condition}_rep{rep1}_vs_rep{rep2}_idr.png"))
            
            # Consensus peaks per condition
            ALL_CONSENSUS_PEAKS.append(os.path.join("consensus_peaks", f"{condition}_consensus.bed"))
            ALL_CONSENSUS_PEAKS.append(os.path.join("consensus_peaks", f"{condition}_consensus_frip.txt"))

# HOMER motif analysis outputs (conditional on run_homer flag)
ALL_HOMER = []
if RUN_HOMER and ALL_CONSENSUS_PEAKS:
    # Extract conditions from consensus peaks
    conditions = set()
    if MACS3_SAMPLES_CONFIG and config.get("metadata_tsv"):
        for sample, info in MACS3_SAMPLES_CONFIG.items():
            if info.get("mark_type"):
                conditions.add(info.get("condition", "default"))
    
    for condition in conditions:
        ALL_HOMER.append(os.path.join("homer_motifs", condition, "homerResults.html"))

# Averaged BigWigs per condition/assay (conditional on metadata TSV)
ALL_AVERAGED_BIGWIGS = []
if MACS3_SAMPLES_CONFIG and config.get("metadata_tsv"):
    # Group by (condition, assay) for averaging
    condition_assay_groups = {}
    for sample, info in MACS3_SAMPLES_CONFIG.items():
        cond = info.get("condition", "default")
        assay = info.get("assay", "unknown")
        key = f"{cond}_{assay}"
        
        if key not in condition_assay_groups:
            condition_assay_groups[key] = []
        condition_assay_groups[key].append(sample)
    
    # Only create averaged BigWigs for groups with ≥2 replicates
    for group_name, samples in condition_assay_groups.items():
        if len(samples) >= 2:
            ALL_AVERAGED_BIGWIGS.append(os.path.join("bigwig_averaged", f"{group_name}_mean.bw"))

# Normalized BigWigs (ChIP/Input log2 ratio) per condition
ALL_NORMALIZED_BIGWIGS = []
if MACS3_SAMPLES_CONFIG and config.get("metadata_tsv"):
    # Get ChIP marks per condition
    chip_conditions = {}
    for sample, info in MACS3_SAMPLES_CONFIG.items():
        if info.get("mark_type"):  # Only ChIP samples
            cond = info.get("condition", "default")
            assay = info.get("assay")
            
            if cond not in chip_conditions:
                chip_conditions[cond] = set()
            chip_conditions[cond].add(assay)
    
    # Create normalized BigWigs for each (condition, chip_mark) pair
    for condition, assays in chip_conditions.items():
        for assay in assays:
            ALL_NORMALIZED_BIGWIGS.append(os.path.join("bigwig_normalized", f"{condition}_{assay}_log2_over_input.bw"))

# Profile plots for averaged/normalized BigWigs per condition
ALL_CONDITION_PROFILES = []
if MACS3_SAMPLES_CONFIG and config.get("metadata_tsv"):
    for condition, assays in chip_conditions.items():
        for assay in assays:
            # Normalized BigWig profiles (log2 ratio)
            ALL_CONDITION_PROFILES.append(os.path.join("condition_profiles", f"{condition}_{assay}_TSS_log2_profile.pdf"))
            ALL_CONDITION_PROFILES.append(os.path.join("condition_profiles", f"{condition}_{assay}_genes_log2_profile.pdf"))
            # Averaged BigWig profiles (mean signal)
            ALL_CONDITION_PROFILES.append(os.path.join("condition_profiles", f"{condition}_{assay}_TSS_mean_profile.pdf"))
            ALL_CONDITION_PROFILES.append(os.path.join("condition_profiles", f"{condition}_{assay}_genes_mean_profile.pdf"))

# =====================================================
# Main Rule
# =====================================================
rule all:
    input:
        ALL_BIGWIGS,
        ALL_FLAGS,
        ALL_SPP,
        MULTIQC_ALIGNMENT,
        ALL_DEEPTOOLS_QC,
        ALL_PROFILE_OUTPUTS,
        ALL_MACS3_PEAKS,
        ALL_FRIP,
        ALL_IDR,
        ALL_CONSENSUS_PEAKS,
        ALL_HOMER,
        ALL_AVERAGED_BIGWIGS,
        ALL_NORMALIZED_BIGWIGS,
        ALL_CONDITION_PROFILES

# =====================================================
# 1. Quality Filtering and Adapter Trimming (Trim_galore)
# =====================================================
rule trim_galore:
    input:
        unpack(get_trim_input)
    output:
        r1=os.path.join(TRIM_DIR, "{sample}_R1_val_1.fq.gz") if SEQ_TYPE == "paired" else os.path.join(TRIM_DIR, "{sample}_trimmed.fq.gz"),
        r2=os.path.join(TRIM_DIR, "{sample}_R2_val_2.fq.gz") if SEQ_TYPE == "paired" else [],
        report1=os.path.join(TRIM_DIR, "{sample}_R1.fastq.gz_trimming_report.txt") if SEQ_TYPE == "paired" else os.path.join(TRIM_DIR, "{sample}.fastq.gz_trimming_report.txt"),
        report2=os.path.join(TRIM_DIR, "{sample}_R2.fastq.gz_trimming_report.txt") if SEQ_TYPE == "paired" else []
    params:
        quality=config["trim_galore"].get("quality", 20),
        min_length=config["trim_galore"].get("min_length", 20),
        adapter=config["trim_galore"].get("adapter", "auto"),
        outdir=TRIM_DIR,
        paired_flag="--paired" if SEQ_TYPE == "paired" else "",
        input_files=lambda wc, input: f"{input.r1} {input.r2}" if SEQ_TYPE == "paired" else input.r1
    threads: 4
    resources:
        mem_mb=config["resources"]["default_mem_mb"]
    log:
        os.path.join(LOG_DIR, "{sample}_trim_galore.log")
    shell:
        """
        trim_galore \
            {params.paired_flag} \
            --{params.adapter} \
            -q {params.quality} \
            --length {params.min_length} \
            --fastqc \
            --cores {threads} \
            -o {params.outdir} \
            {params.input_files} \
            > {log} 2>&1
        """

# =====================================================
# 2. Automatic Genome File Download (Optional)
# =====================================================
# Download and decompress genome FASTA if URL provided
if config["genome"].get("fasta_url"):
    rule download_genome_fasta:
        output:
            temp(GENOME_FASTA + ".gz")
        params:
            url=config["genome"]["fasta_url"]
        log:
            os.path.join(LOG_DIR, "download_genome_fasta.log")
        shell:
            """
            echo "Downloading genome FASTA from {params.url}..."
            wget -O {output} {params.url} > {log} 2>&1
            """
    
    rule decompress_fasta:
        input:
            GENOME_FASTA + ".gz"
        output:
            GENOME_FASTA
        log:
            os.path.join(LOG_DIR, "decompress_fasta.log")
        shell:
            """
            echo "Decompressing genome FASTA..."
            gunzip -c {input} > {output} 2> {log}
            """

# Download and decompress GTF if URL provided
if config["genome"].get("gtf_url"):
    rule download_gtf:
        output:
            temp(GTF_FILE + ".gz")
        params:
            url=config["genome"]["gtf_url"]
        log:
            os.path.join(LOG_DIR, "download_gtf.log")
        shell:
            """
            echo "Downloading GTF from {params.url}..."
            wget -O {output} {params.url} > {log} 2>&1
            """
    
    rule decompress_gtf:
        input:
            GTF_FILE + ".gz"
        output:
            GTF_FILE
        log:
            os.path.join(LOG_DIR, "decompress_gtf.log")
        shell:
            """
            echo "Decompressing GTF..."
            gunzip -c {input} > {output} 2> {log}
            """

# =====================================================
# 3. Genome Indexing (Bowtie2)
# =====================================================
rule build_index:
    input:
        GENOME_FASTA
    output:
        expand(
            os.path.join(os.path.dirname(GENOME_INDEX) or ".", "{{prefix}}.{ext}.bt2"),
            ext=["1", "2", "3", "4", "rev.1", "rev.2"]
        )
    params:
        prefix=GENOME_INDEX
    threads: THREADS
    resources:
        mem_mb=config["resources"]["default_mem_mb"]
    log:
        os.path.join(LOG_DIR, "bowtie2_build.log")
    shell:
        """
        mkdir -p $(dirname {params.prefix})
        bowtie2-build --threads {threads} {input} {params.prefix} > {log} 2>&1
        """

# =====================================================
# 3. Alignment, Filtering, and Sorting
# =====================================================
rule align_sort:
    input:
        unpack(get_align_input),
        index=expand(
            os.path.join(os.path.dirname(GENOME_INDEX) or ".", os.path.basename(GENOME_INDEX) + ".{ext}.bt2"),
            ext=["1", "2", "3", "4", "rev.1", "rev.2"]
        )
    output:
        bam=os.path.join(OUT_ALIGN, "{sample}.sorted.bam"),
        bai=os.path.join(OUT_ALIGN, "{sample}.sorted.bam.bai")
    params:
        index=GENOME_INDEX,
        # Bowtie2 parameters from config
        sensitivity=config["bowtie2"].get("sensitivity", "--very-sensitive"),
        no_unal="--no-unal" if config["bowtie2"].get("no_unal", True) else "",
        no_mixed="--no-mixed" if SEQ_TYPE == "paired" and config["bowtie2"].get("no_mixed", True) else "",
        no_discordant="--no-discordant" if SEQ_TYPE == "paired" and config["bowtie2"].get("no_discordant", True) else "",
        # Samtools filter parameters from config
        mapq=config["samtools"].get("mapq_threshold", 30),
        flags=config["samtools"].get("sam_flags_exclude", 768)
    threads: THREADS
    resources:
        mem_mb=config["resources"]["align_mem_mb"]
    log:
        bowtie2=os.path.join(LOG_DIR, "{sample}_bowtie2.log"),
        samtools_sort=os.path.join(LOG_DIR, "{sample}_samtools_sort.log")
    run:
        if SEQ_TYPE == "paired":
            # Paired-end alignment
            shell("""
                bowtie2 \
                  -p {threads} \
                  -x {params.index} \
                  -1 {input.r1} \
                  -2 {input.r2} \
                  {params.sensitivity} \
                  {params.no_mixed} \
                  {params.no_discordant} \
                  {params.no_unal} \
                  2> {log.bowtie2} \
                | samtools view -@ {threads} -bS -q {params.mapq} -F {params.flags} - \
                | samtools sort -@ {threads} -o {output.bam} - 2> {log.samtools_sort}
                
                samtools index -@ {threads} {output.bam}
            """)
        else:
            # Single-end alignment
            shell("""
                bowtie2 \
                  -p {threads} \
                  -x {params.index} \
                  -U {input.r1} \
                  {params.sensitivity} \
                  {params.no_unal} \
                  2> {log.bowtie2} \
                | samtools view -@ {threads} -bS -q {params.mapq} -F {params.flags} - \
                | samtools sort -@ {threads} -o {output.bam} - 2> {log.samtools_sort}
                
                samtools index -@ {threads} {output.bam}
            """)

# Interactive FRiP table with Datavzrd
rule view_frip_with_datavzrd:
    input:
        config="datavzrd/frip_summary.yaml",
        table="frip_analysis/all_samples_frip.tsv"
    output:
        report(
            directory("frip_analysis/interactive_table"),
            htmlindex="index.html",
            caption="report/frip_summary.rst",
            category="Quality Control",
            labels={"type": "FRiP summary"}
        )
    log:
        "logs/datavzrd_frip.log"
    wrapper:
        "v4.7.2/utils/datavzrd"

# =====================================================
# 4. Duplicate Removal (Proper Fixmate Workflow)
# =====================================================
rule markdup:
    input:
        bam=os.path.join(OUT_ALIGN, "{sample}.sorted.bam")
    output:
        bam=os.path.join(OUT_DEDUP, "{sample}.dedup.bam")
    resources:
        mem_mb=config["resources"]["default_mem_mb"]
    log:
        os.path.join(LOG_DIR, "{sample}_markdup.log")
    threads: 1
    shell:
        """
        # Unified 4-step fixmate workflow (works for both SE and PE)
        # Step 1: Collate (group read pairs / reads by name)
        # Step 2: Fixmate (add mate tags - only useful for PE, harmless for SE)
        # Step 3: Sort by coordinate
        # Step 4: Mark duplicates
        
        samtools collate -@ {threads} -u {input.bam} | \\
          samtools fixmate -m -u stdin stdout | \\
          samtools sort -@ {threads} -u stdin | \\
          samtools markdup -@ {threads} stdin {output.bam} \\
          2> {log}
        
        samtools index {output.bam}
        """

# =====================================================
# 5. Quality Control Metrics
# =====================================================
rule qc_flagstat_stats_pbc:
    input:
        bam=os.path.join(OUT_DEDUP, "{sample}.dedup.bam")
    output:
        flagstat=os.path.join(OUT_QC, "{sample}.flagstat.txt"),
        stats=os.path.join(OUT_QC, "{sample}.stats.txt"),
        read5=temp(os.path.join(OUT_QC, "{sample}.read5.bed")),
        pbc=os.path.join(OUT_QC, "{sample}.pbc.txt")
    threads: 1
    resources:
        mem_mb=config["resources"]["default_mem_mb"]
    log:
        os.path.join(LOG_DIR, "{sample}_qc.log")
    shell:
        """
        # Flagstat
        samtools flagstat {input.bam} > {output.flagstat} 2> {log}
        
        # Stats
        samtools stats {input.bam} > {output.stats} 2>> {log}
        
        # PBC metrics
        bedtools bamtobed -i {input.bam} \
        | awk 'BEGIN{{OFS="\\t"}} ($6=="+"){{print $1,$2,$2+1}} ($6=="-"){{print $1,$3-1,$3}}' \
        | sort -k1,1 -k2,2n > {output.read5}
        
        uniq -c {output.read5} \
        | awk '{{c=$1; total+=c; uniq++; if(c==1) single++; if(c==2) double++;}} 
               END{{
                   if(total>0){{
                       NRF=uniq/total;
                       PBC1=single/uniq;
                       if(double) PBC2=single/double; else PBC2="Inf";
                   }} else {{
                       NRF=0; PBC1=0; PBC2=0
                   }}
                   printf("NRF=%.3f\\tPBC1=%.3f\\tPBC2=%s\\n", NRF, PBC1, PBC2);
               }}' > {output.pbc}
        """

# =====================================================
# 6. MultiQC - Alignment QC Report
# =====================================================
rule multiqc_alignment:
    input:
        flagstats=expand(os.path.join(OUT_QC, "{sample}.flagstat.txt"), sample=SAMPLES),
        stats=expand(os.path.join(OUT_QC, "{sample}.stats.txt"), sample=SAMPLES)
    output:
        report(
            os.path.join(OUT_QC, "alignment_qc_report.html"),
            caption="report/multiqc_alignment.rst",
            category="Quality Control",
            labels={"type": "alignment metrics"}
        ),
        data=directory(os.path.join(OUT_QC, "alignment_qc_report_data"))
    params:
        outdir=OUT_QC,
        name="alignment_qc_report"
    log:
        os.path.join(LOG_DIR, "multiqc_alignment.log")
    shell:
        """
        multiqc {params.outdir} \
            -o {params.outdir} \
            -n {params.name} \
            --force \
            > {log} 2>&1
        """

# =====================================================
# 7. Optional SPP Analysis
# =====================================================
if RUN_SPP:
    rule spp_qc:
        input:
            bam=os.path.join(OUT_DEDUP, "{sample}.dedup.bam")
        output:
            avp=os.path.join(OUT_QC, "{sample}_avp.pdf"),
            sppqc=os.path.join(OUT_QC, "{sample}_spp.qc.txt")
        threads: 1
        resources:
            mem_mb=config["resources"]["default_mem_mb"]
        log:
            os.path.join(LOG_DIR, "{sample}_spp.log")
        shell:
            """
            Rscript {SPP_SCRIPT} \
                -c={input.bam} \
                -savp={output.avp} \
                -out={output.sppqc} \
                > {log} 2>&1 || echo "SPP failed for {wildcards.sample}" >> {log}
            """

# =====================================================
# 8. DeepTools QC (Multi-sample)
# =====================================================
rule multiBamSummary:
    input:
        bams=expand(os.path.join(OUT_DEDUP, "{sample}.dedup.bam"), sample=SAMPLES)
    output:
        npz=os.path.join(OUT_DEEPTOOLS, "all_samples_bins.npz"),
        raw=os.path.join(OUT_DEEPTOOLS, "all_samples_bins.tab")
    params:
        labels=" ".join(SAMPLES)
    threads: THREADS
    resources:
        mem_mb=config["resources"]["deeptools_mem_mb"]
    log:
        os.path.join(LOG_DIR, "multiBamSummary.log")
    shell:
        """
        multiBamSummary bins \
            -b {input.bams} \
            --labels {params.labels} \
            -o {output.npz} \
            --outRawCounts {output.raw} \
            --numberOfProcessors {threads} \
            > {log} 2>&1
        """

rule plotCorrelation:
    input:
        npz=os.path.join(OUT_DEEPTOOLS, "all_samples_bins.npz")
    output:
        heatmap=report(
            os.path.join(OUT_DEEPTOOLS, "correlation_heatmap.pdf"),
            caption="report/correlation_plot.rst",
            category="QC Plots",
            labels={"plot": "correlation_heatmap"}
        ),
        scatter=report(
            os.path.join(OUT_DEEPTOOLS, "correlation_scatterplot.pdf"),
            caption="report/correlation_plot.rst",
            category="QC Plots",
            labels={"plot": "correlation_scatter"}
        ),
        matrix=os.path.join(OUT_DEEPTOOLS, "correlation_matrix.tab")
    params:
        dpi=config.get("deeptools_plots", {}).get("dpi", 200),
        width=config.get("deeptools_plots", {}).get("correlation", {}).get("plot_width", 11),
        height=config.get("deeptools_plots", {}).get("correlation", {}).get("plot_height", 9.5)
    log:
        os.path.join(LOG_DIR, "plotCorrelation.log")
    shell:
        """
        plotCorrelation \
            -in {input.npz} \
            --corMethod spearman \
            --skipZeros \
            --whatToPlot heatmap \
            --colorMap RdYlBu \
            --plotNumbers \
            --plotTitle "Spearman Correlation" \
            --plotWidth {params.width} \
            --plotHeight {params.height} \
            --dpi {params.dpi} \
            -o {output.heatmap} \
            --outFileCorMatrix {output.matrix} \
            > {log} 2>&1
        """

rule plotPCA:
    input:
        npz=os.path.join(OUT_DEEPTOOLS, "all_samples_bins.npz")
    output:
        report(
            os.path.join(OUT_DEEPTOOLS, "pca.pdf"),
            caption="report/pca_plot.rst",
            category="QC Plots",
            labels={"plot": "PCA"}
        ),
        tab=os.path.join(OUT_DEEPTOOLS, "PCA.tab")
    params:
        dpi=config.get("deeptools_plots", {}).get("dpi", 200),
        width=config.get("deeptools_plots", {}).get("pca", {}).get("plot_width", 14),
        height=config.get("deeptools_plots", {}).get("pca", {}).get("plot_height", 12),
        markers=config.get("deeptools_plots", {}).get("pca", {}).get("markers", []),
        colors=config.get("deeptools_plots", {}).get("pca", {}).get("colors", [])
    log:
        os.path.join(LOG_DIR, "plotPCA.log")
    run:
        # Build command with optional markers and colors
        cmd = f"""
        plotPCA \\
            -in {input.npz} \\
            -o {output[0]} \\
            -T "PCA of Binned Coverage" \\
            --transpose \\
            --plotWidth {params.width} \\
            --plotHeight {params.height} \\
            --dpi {params.dpi} \\
            --outFileNameData {output.tab}"""
        
        # Add markers if specified
        if params.markers:
            markers_str = " ".join([f"'{m}'" for m in params.markers])
            cmd += f" \\\n            --markers {markers_str}"
        
        # Add colors if specified
        if params.colors:
            colors_str = " ".join([f"'{c}'" for c in params.colors])
            cmd += f" \\\n            --colors {colors_str}"
        
        cmd += f" \\\n            > {log} 2>&1"
        
        shell(cmd)

rule plotFingerprint:
    input:
        bams=expand(os.path.join(OUT_DEDUP, "{sample}.dedup.bam"), sample=SAMPLES)
    output:
        fingerprints=report(
            os.path.join(OUT_DEEPTOOLS, "fingerprints.pdf"),
            caption="report/pca_plot.rst",
            category="QC Plots",
            labels={"plot": "fingerprint"}
        )
    params:
        labels=" ".join(SAMPLES)
    threads: THREADS
    resources:
        mem_mb=config["resources"]["deeptools_mem_mb"]
    log:
        os.path.join(LOG_DIR, "plotFingerprint.log")
    shell:
        """
        plotFingerprint \
            -b {input.bams} \
            --labels {params.labels} \
            --skipZeros \
            --numberOfSamples 50000 \
            -T "Fingerprints" \
            --plotFile {output.pdf} \
            --numberOfProcessors {threads} \
            > {log} 2>&1
        """

rule plotCoverage:
    input:
        bams=expand(os.path.join(OUT_DEDUP, "{sample}.dedup.bam"), sample=SAMPLES)
    output:
        pdf=os.path.join(OUT_DEEPTOOLS, "coverage_histogram.pdf"),
        counts=os.path.join(OUT_DEEPTOOLS, "coverage_counts.txt")
    params:
        labels=" ".join(SAMPLES)
    threads: THREADS
    resources:
        mem_mb=config["resources"]["deeptools_mem_mb"]
    log:
        os.path.join(LOG_DIR, "plotCoverage.log")
    shell:
        """
        plotCoverage \
            -b {input.bams} \
            --labels {params.labels} \
            --numberOfSamples 1000000 \
            --ignoreDuplicates \
            --minMappingQuality 30 \
            -o {output.pdf} \
            --outRawCounts {output.counts} \
            --numberOfProcessors {threads} \
            > {log} 2>&1
        """

# =====================================================
# 10. BigWig Coverage Tracks
# =====================================================
rule bamCoverage:
    input:
        bam=os.path.join(OUT_DEDUP, "{sample}.dedup.bam")
    output:
        bw=os.path.join(OUT_DEEPTOOLS, "{sample}.dedup.bw")
    params:
        normalization=f"--normalizeUsing {NORMALIZATION_METHOD}" if NORMALIZATION_METHOD and NORMALIZATION_METHOD != "None" else ""
    threads: THREADS
    resources:
        mem_mb=config["resources"]["default_mem_mb"]
    log:
        os.path.join(LOG_DIR, "{sample}_bamCoverage.log")
    shell:
        """
        bamCoverage \
            -b {input.bam} \
            -o {output.bw} \
            {params.normalization} \
            --effectiveGenomeSize {EFFECTIVE_GENOME_SIZE} \
            --binSize {BIN_SIZE} \
            --smoothLength {SMOOTH_LENGTH} \
            --numberOfProcessors {threads} \
            --ignoreDuplicates \
            > {log} 2>&1
        """

# =====================================================
# 11. Prepare Genomic Regions from GTF
# =====================================================
rule prepare_regions:
    input:
        gtf=GTF_FILE
    output:
        tss=os.path.join(OUT_DEEPTOOLS, "tss.bed"),
        genes=os.path.join(OUT_DEEPTOOLS, "genes.bed")
    threads: 1
    log:
        os.path.join(LOG_DIR, "prepare_regions.log")
    shell:
        """
        # Extract TSS positions
        awk 'BEGIN{{OFS="\\t"}} $3=="transcript" {{
            if($7=="+") print $1, $4-1, $4, $12, ".", $7;
            else print $1, $5-1, $5, $12, ".", $7
        }}' {input.gtf} \
        | tr -d '";' \
        | sort -k1,1V -k2,2n > {output.tss} 2> {log}
        
        # Extract gene bodies
        awk 'BEGIN{{OFS="\\t"}} $3=="gene" {{
            print $1, $4-1, $5, $10, ".", $7
        }}' {input.gtf} \
        | tr -d '";' \
        | sort -k1,1V -k2,2n > {output.genes} 2>> {log}
        """

# =====================================================
# 12. Profile Plots (TSS and Gene Body)
# =====================================================
rule computeMatrix_TSS:
    input:
        tss=os.path.join(OUT_DEEPTOOLS, "tss.bed"),
        bws=expand(os.path.join(OUT_DEEPTOOLS, "{sample}.dedup.bw"), sample=SAMPLES)
    output:
        matrix=os.path.join(OUT_DEEPTOOLS, "matrix_all_bw_TSS.gz")
    params:
        upstream=config["parameters"]["tss_upstream"],
        downstream=config["parameters"]["tss_downstream"]
    threads: THREADS
    resources:
        mem_mb=config["resources"]["deeptools_mem_mb"]
    log:
        os.path.join(LOG_DIR, "computeMatrix_TSS.log")
    shell:
        """
        computeMatrix reference-point \
            --referencePoint TSS \
            -b {params.upstream} \
            -a {params.downstream} \
            -R {input.tss} \
            -S {input.bws} \
            --skipZeros \
            -o {output.matrix} \
            --binSize {BIN_SIZE} \
            --numberOfProcessors {threads} \
            > {log} 2>&1
        """

rule plotProfile_TSS:
    input:
        matrix=os.path.join(OUT_DEEPTOOLS, "matrix_all_bw_TSS.gz")
    output:
        profile=os.path.join(OUT_DEEPTOOLS, "all_samples_TSS_profile.pdf")
    log:
        os.path.join(LOG_DIR, "plotProfile_TSS.log")
    shell:
        """
        plotProfile \
            --matrixFile {input.matrix} \
            --outFileName {output.profile} \
            --perGroup \
            --refPointLabel TSS \
            --plotType se \
            --legendLocation upper-right \
            > {log} 2>&1
        """

rule computeMatrix_body:
    input:
        genes=os.path.join(OUT_DEEPTOOLS, "genes.bed"),
        bws=expand(os.path.join(OUT_DEEPTOOLS, "{sample}.dedup.bw"), sample=SAMPLES)
    output:
        matrix=os.path.join(OUT_DEEPTOOLS, "matrix_all_bw_scalar.gz")
    params:
        body_length=config["parameters"]["gene_body_length"],
        upstream=config["parameters"]["gene_body_upstream"],
        downstream=config["parameters"]["gene_body_downstream"]
    threads: THREADS
    resources:
        mem_mb=config["resources"]["deeptools_mem_mb"]
    log:
        os.path.join(LOG_DIR, "computeMatrix_body.log")
    shell:
        """
        computeMatrix scale-regions \
            -R {input.genes} \
            -S {input.bws} \
            --regionBodyLength {params.body_length} \
            -b {params.upstream} \
            -a {params.downstream} \
            --binSize {BIN_SIZE} \
            --skipZeros \
            -o {output.matrix} \
            --numberOfProcessors {threads} \
            > {log} 2>&1
        """

rule plotProfile_body:
    input:
        matrix=os.path.join(OUT_DEEPTOOLS, "matrix_all_bw_scalar.gz")
    output:
        profile=os.path.join(OUT_DEEPTOOLS, "all_samples_scalar_profile.pdf")
    log:
        os.path.join(LOG_DIR, "plotProfile_body.log")
    shell:
        """
        plotProfile \
            --matrixFile {input.matrix} \
            --outFileName {output.profile} \
            --perGroup \
            --plotType se \
            --legendLocation upper-right \
            > {log} 2>&1
        """

# =====================================================
# 13. MACS3 Peak Calling (Per-Replicate)
# =====================================================
if MACS3_SAMPLES:
    rule macs3_callpeak:
        input:
            treatment=os.path.join(OUT_DEDUP, "{sample}.dedup.bam"),
            control=lambda wc: os.path.join(OUT_DEDUP, f"{MACS3_SAMPLES_CONFIG[wc.sample]['control']}.dedup.bam") 
                    if MACS3_SAMPLES_CONFIG[wc.sample].get("control") else []
        output:
            xls=os.path.join("macs3_results", "{sample}_peaks.xls"),
            narrowPeak=os.path.join("macs3_results", "{sample}_peaks.narrowPeak"),
            summits=os.path.join("macs3_results", "{sample}_summits.bed"),
            broadPeak=os.path.join("macs3_results", "{sample}_peaks.broadPeak"),
            gappedPeak=os.path.join("macs3_results", "{sample}_peaks.gappedPeak")
        params:
            name="{sample}",
            outdir="macs3_results",
            genome_size=config.get("macs3_params", {}).get("genome_size", "hs"),
            pvalue=config.get("macs3_params", {}).get("pvalue_cutoff", 0.01),
            keep_dup=config.get("macs3_params", {}).get("keep_dup", "all"),
            broad_cutoff=config.get("macs3_params", {}).get("broad_cutoff", 0.1),
            mark_type=lambda wc: MACS3_SAMPLES_CONFIG[wc.sample]["mark_type"],
            control_flag=lambda wc, input: f"-c {input.control}" if input.control else "",
            # Optional parameters
            format=config.get("macs3_params", {}).get("format", "BAM"),
            bdg_flag="-B" if config.get("macs3_params", {}).get("generate_bedgraph", False) else "",
            spmr_flag="--SPMR" if config.get("macs3_params", {}).get("bedgraph_spmr", False) else "",
            model_flag="--nomodel --extsize {}".format(config.get("macs3_params", {}).get("extsize", 200)) 
                       if config.get("macs3_params", {}).get("no_model", False) else ""
        log:
            os.path.join(LOG_DIR, "{sample}_macs3.log")
        wildcard_constraints:
            sample="|".join(MACS3_SAMPLES)
        run:
            # Build MACS3 command based on mark type
            broad_flag = f"--broad --broad-cutoff {params.broad_cutoff}" if params.mark_type == "broad" else ""
            
            shell(f"""
                mkdir -p {params.outdir}
                
                macs3 callpeak \\
                    -t {input.treatment} \\
                    {params.control_flag} \\
                    -f {params.format} \\
                    -g {params.genome_size} \\
                    -n {params.name} \\
                    -p {params.pvalue} \\
                    --keep-dup {params.keep_dup} \\
                    {broad_flag} \\
                    {params.bdg_flag} \\
                    {params.spmr_flag} \\
                    {params.model_flag} \\
                    --outdir {params.outdir} \\
                    > {log} 2>&1
            """)

# =====================================================
# 14. FRiP (Fraction of Reads in Peaks) Calculation
# =====================================================
if MACS3_SAMPLES:
    rule calculate_frip:
        input:
            bam=os.path.join(OUT_DEDUP, "{sample}.dedup.bam"),
            peaks=lambda wc: os.path.join("macs3_results", f"{wc.sample}_peaks.narrowPeak") 
                   if MACS3_SAMPLES_CONFIG[wc.sample]["mark_type"] == "narrow"
                   else os.path.join("macs3_results", f"{wc.sample}_peaks.broadPeak")
        output:
            frip=os.path.join("frip_analysis", "{sample}.frip.txt"),
            merged_peaks=os.path.join("frip_analysis", "{sample}.peaks.merged.bed")
        log:
            os.path.join(LOG_DIR, "{sample}_frip.log")
        wildcard_constraints:
            sample="|".join(MACS3_SAMPLES)
        shell:
            """
            mkdir -p frip_analysis
            
            # Count total mapped reads
            TOTAL_READS=$(samtools view -c {input.bam})
            
            # Merge overlapping peaks to avoid double-counting
            sort -k1,1 -k2,2n {input.peaks} | \\
              bedtools merge -i stdin > {output.merged_peaks}
            
            # Count reads overlapping peaks
            READS_IN_PEAKS=$(samtools view -b {input.bam} | \\
              bedtools intersect -u -a stdin -b {output.merged_peaks} | \\
              samtools view -c)
            
            # Calculate FRiP
            FRIP=$(awk -v n="$READS_IN_PEAKS" -v d="$TOTAL_READS" 'BEGIN {{printf "%.5f", n/d}}')
            
            # Write results
            echo "Sample: {wildcards.sample}" > {output.frip}
            echo "Total mapped reads : $TOTAL_READS" >> {output.frip}
            echo "Reads in peaks     : $READS_IN_PEAKS" >> {output.frip}
            echo "FRiP               : $FRIP" >> {output.frip}
            
            echo "FRiP calculation completed for {wildcards.sample}" > {log} 2>&1
            """
    
    rule aggregate_frip:
        input:
            frip_files=expand(os.path.join("frip_analysis", "{sample}.frip.txt"), sample=MACS3_SAMPLES)
        output:
            summary=os.path.join("frip_analysis", "all_samples_frip.tsv")
        log:
            os.path.join(LOG_DIR, "aggregate_frip.log")
        run:
            import re
            
            with open(output.summary, 'w') as out:
                out.write("Sample\tTotal_Reads\tReads_in_Peaks\tFRiP\n")
                
                for frip_file in input.frip_files:
                    with open(frip_file) as f:
                        content = f.read()
                        sample = re.search(r'Sample: (\S+)', content).group(1)
                        total = re.search(r'Total mapped reads\s+:\s+(\d+)', content).group(1)
                        in_peaks = re.search(r'Reads in peaks\s+:\s+(\d+)', content).group(1)
                        frip = re.search(r'FRiP\s+:\s+([\d.]+)', content).group(1)
                        
                        out.write(f"{sample}\t{total}\t{in_peaks}\t{frip}\n")
            
            shell("echo 'FRiP summary table created' > {log} 2>&1")

# =====================================================
# 15. IDR Analysis (Irreproducible Discovery Rate)
# =====================================================
if MACS3_SAMPLES_CONFIG and config.get("metadata_tsv"):
    rule idr_pairwise:
        input:
            peak1=lambda wc: os.path.join("macs3_results", f"{wc.sample1}_peaks.{wc.peak_type}Peak"),
            peak2=lambda wc: os.path.join("macs3_results", f"{wc.sample2}_peaks.{wc.peak_type}Peak")
        output:
            idr_output=os.path.join("idr_analysis", "{condition}_rep{rep1}_vs_rep{rep2}_idr.txt"),
            idr_plot=os.path.join("idr_analysis", "{condition}_rep{rep1}_vs_rep{rep2}_idr.png")
        params:
            peak_type="{peak_type}"
        log:
            os.path.join(LOG_DIR, "{condition}_rep{rep1}_vs_rep{rep2}_idr.log")
        shell:
            """
            mkdir -p idr_analysis
            
            idr --samples {input.peak1} {input.peak2} \\
                --input-file-type {params.peak_type}Peak \\
                --rank signal.value \\
                --output-file {output.idr_output} \\
                --plot \\
                --log-output-file {log}
            """
    
    rule consensus_peaks:
        input:
            idr_files=lambda wc: expand(
                os.path.join("idr_analysis", f"{{condition}}_rep{{rep1}}_vs_rep{{rep2}}_idr.txt"),
                condition=[wc.condition],
                rep1=[r1 for r1, r2 in get_replicate_pairs(wc.condition)],
                rep2=[r2 for r1, r2 in get_replicate_pairs(wc.condition)]
            ),
            bam=lambda wc: os.path.join(OUT_DEDUP, f"{get_first_sample(wc.condition)}.dedup.bam")
        output:
            bed=report(
                os.path.join("consensus_peaks", "{condition}_consensus.bed"),
                caption="report/consensus_peaks.rst",
                category="Consensus Peaks",
                labels={"condition": "{condition}"}
            ),
            frip=os.path.join("consensus_peaks", "{condition}_consensus_frip.txt")
        log:
            os.path.join(LOG_DIR, "{condition}_consensus.log")
        shell:
            """
            mkdir -p consensus_peaks
            
            # Filter peaks with IDR >= 1.3 (more stringent threshold)
            # Merge all passing peaks from all pairwise comparisons
            cat {input.idr_files} | awk '$12 >= 1.3 {{print $1"\\t"$2"\\t"$3}}' | \\
              sort -k1,1 -k2,2n | \\
              bedtools merge -i stdin > {output.consensus}
            
            # Calculate FRiP on consensus peaks
            TOTAL_READS=$(samtools view -c {input.bam})
            READS_IN_PEAKS=$(samtools view -b {input.bam} | \\
              bedtools intersect -u -a stdin -b {output.consensus} | \\
              samtools view -c)
            FRIP=$(awk -v n="$READS_IN_PEAKS" -v d="$TOTAL_READS" 'BEGIN {{printf "%.5f", n/d}}')
            
            echo "Condition: {wildcards.condition}" > {output.frip}
            echo "Total reads: $TOTAL_READS" >> {output.frip}
            echo "Reads in consensus peaks: $READS_IN_PEAKS" >> {output.frip}
            echo "FRiP (consensus): $FRIP" >> {output.frip}
            echo "Number of consensus peaks: $(wc -l < {output.consensus})" >> {output.frip}
            
            echo "Consensus peaks generated for {wildcards.condition}" > {log} 2>&1
            """

# Helper functions for IDR
def get_replicate_pairs(condition):
    """Get all pairwise replicate combinations for a condition"""
    from itertools import combinations
    
    samples = [s for s, info in MACS3_SAMPLES_CONFIG.items() 
               if info.get("mark_type") and info.get("condition") == condition]
    
    reps = sorted([MACS3_SAMPLES_CONFIG[s]["replicate"] for s in samples])
    return list(combinations(reps, 2))

def get_first_sample(condition):
    """Get first ChIP sample name for a condition (for FRiP calculation)"""
    for s, info in MACS3_SAMPLES_CONFIG.items():
        if info.get("mark_type") and info.get("condition") == condition:
            return s
    return None

# =====================================================
# 16. HOMER Motif Analysis (Optional)
# =====================================================
if RUN_HOMER and MACS3_SAMPLES_CONFIG and config.get("metadata_tsv"):
    rule homer_motif_analysis:
        input:
            consensus=os.path.join("consensus_peaks", "{condition}_consensus.bed"),
            genome=GENOME_FASTA
        output:
            html=os.path.join("homer_motifs", "{condition}", "homerResults.html"),
            motifs=directory(os.path.join("homer_motifs", "{condition}"))
        params:
            size=config.get("homer_params", {}).get("size", 200),
            length=config.get("homer_params", {}).get("length", "8,10,12"),
            mask="-mask" if config.get("homer_params", {}).get("mask", True) else ""
        threads: THREADS
        log:
            os.path.join(LOG_DIR, "{condition}_homer.log")
        shell:
            """
            mkdir -p {output.motifs}
            
            findMotifsGenome.pl \
                {input.consensus} \
                {input.genome} \
                {output.motifs} \
                -size {params.size} \
                -len {params.length} \
                {params.mask} \
                -p {threads} \
                > {log} 2>&1
            """

# =====================================================
# 17. Average BigWigs Per Condition/Assay
# =====================================================
if MACS3_SAMPLES_CONFIG and config.get("metadata_tsv"):
    rule average_bigwigs:
        input:
            bigwigs=lambda wc: expand(
                os.path.join(OUT_DEEPTOOLS, "{sample}.dedup.bw"),
                sample=get_samples_for_group(wc.condition, wc.assay)
            )
        output:
            mean_bw=os.path.join("bigwig_averaged", "{condition}_{assay}_mean.bw")
        threads: THREADS
        log:
            os.path.join(LOG_DIR, "{condition}_{assay}_bigwigAverage.log")
        shell:
            """
            mkdir -p bigwig_averaged
            
            bigwigAverage \
                -b {input.bigwigs} \
                -o {output.mean_bw} \
                -p {threads} \
                > {log} 2>&1
            """

# Helper function for BigWig averaging
def get_samples_for_group(condition, assay):
    """Get all sample names for a given (condition, assay) combination"""
    samples = []
    for sample, info in MACS3_SAMPLES_CONFIG.items():
        if info.get("condition") == condition and info.get("assay") == assay:
            samples.append(sample)
    return samples

# =====================================================
# 18. Normalize BigWigs (ChIP/Input Log2 Ratio)
# =====================================================
if MACS3_SAMPLES_CONFIG and config.get("metadata_tsv"):
    rule normalize_bigwigs:
        input:
            chip_bw=os.path.join("bigwig_averaged", "{condition}_{assay}_mean.bw"),
            input_bw=os.path.join("bigwig_averaged", "{condition}_input_mean.bw")
        output:
            normalized_bw=os.path.join("bigwig_normalized", "{condition}_{assay}_log2_over_input.bw")
        params:
            pseudocount=1
        threads: THREADS
        log:
            os.path.join(LOG_DIR, "{condition}_{assay}_bigwigCompare.log")
        shell:
            """
            mkdir -p bigwig_normalized
            
            bigwigCompare \
                -b1 {input.chip_bw} \
                -b2 {input.input_bw} \
                --operation log2 \
                --pseudocount {params.pseudocount} \
                -p {threads} \
                -o {output.normalized_bw} \
                > {log} 2>&1
            """

# =====================================================
# 19. Profile Plots for Averaged/Normalized BigWigs
# =====================================================
if MACS3_SAMPLES_CONFIG and config.get("metadata_tsv"):
    # TSS profiles - normalized (log2)
    rule profile_tss_normalized:
        input:
            tss=os.path.join(OUT_DEEPTOOLS, "tss.bed"),
            bw=os.path.join("bigwig_normalized", "{condition}_{assay}_log2_over_input.bw")
        output:
            matrix=os.path.join("condition_profiles", "{condition}_{assay}_TSS_log2.mat.gz"),
            profile=report(
                os.path.join("condition_profiles", "{condition}_{assay}_TSS_log2_profile.pdf"),
                caption="report/tss_profile.rst",
                category="Profile Plots",
                labels={"condition": "{condition}", "assay": "{assay}", "type": "TSS_normalized"}
            ),
            pdf=os.path.join("condition_profiles", "{condition}_{assay}_TSS_log2_profile.pdf")
        threads: THREADS
        log:
            os.path.join(LOG_DIR, "{condition}_{assay}_TSS_log2_profile.log")
        shell:
            """
            mkdir -p condition_profiles
            
            computeMatrix reference-point \
                --referencePoint TSS \
                -b 3000 -a 3000 \
                --binSize 250 \
                -R {input.tss} \
                -S {input.bw} \
                --skipZeros --missingDataAsZero \
                -p {threads} \
                -o {output.matrix} \
                > {log} 2>&1
            
            plotProfile \
                -m {output.matrix} \
                --refPointLabel TSS \
                --yAxisLabel "log2(IP/Input)" \
                --plotTitle "{wildcards.condition} {wildcards.assay} log2(IP/Input)" \
                -out {output.pdf} \
                >> {log} 2>&1
            """
    
    # Gene body profiles - normalized (log2)
    rule profile_genes_normalized:
        input:
            genes=os.path.join(OUT_DEEPTOOLS, "genes.bed"),
            bw=os.path.join("bigwig_normalized", "{condition}_{assay}_log2_over_input.bw")
        output:
            matrix=os.path.join("condition_profiles", "{condition}_{assay}_genes_log2.mat.gz"),
            pdf=os.path.join("condition_profiles", "{condition}_{assay}_genes_log2_profile.pdf")
        threads: THREADS
        log:
            os.path.join(LOG_DIR, "{condition}_{assay}_genes_log2_profile.log")
        shell:
            """
            mkdir -p condition_profiles
            
            computeMatrix scale-regions \
                -b 3000 -a 3000 \
                --regionBodyLength 5000 \
                --binSize 250 \
                -R {input.genes} \
                -S {input.bw} \
                --skipZeros --missingDataAsZero \
                -p {threads} \
                -o {output.matrix} \
                > {log} 2>&1
            
            plotProfile \
                -m {output.matrix} \
                --yAxisLabel "log2(IP/Input)" \
                --plotTitle "{wildcards.condition} {wildcards.assay} gene-body log2(IP/Input)" \
                -out {output.pdf} \
                >> {log} 2>&1
            """
    
    # TSS profiles - averaged (mean signal)
    rule profile_tss_mean:
        input:
            tss=os.path.join(OUT_DEEPTOOLS, "tss.bed"),
            bw=os.path.join("bigwig_averaged", "{condition}_{assay}_mean.bw")
        output:
            matrix=os.path.join("condition_profiles", "{condition}_{assay}_TSS_mean.mat.gz"),
            pdf=os.path.join("condition_profiles", "{condition}_{assay}_TSS_mean_profile.pdf")
        threads: THREADS
        log:
            os.path.join(LOG_DIR, "{condition}_{assay}_TSS_mean_profile.log")
        shell:
            """
            mkdir -p condition_profiles
            
            computeMatrix reference-point \
                --referencePoint TSS \
                -b 3000 -a 3000 \
                --binSize 250 \
                -R {input.tss} \
                -S {input.bw} \
                --skipZeros --missingDataAsZero \
                -p {threads} \
                -o {output.matrix} \
                > {log} 2>&1
            
            plotProfile \
                -m {output.matrix} \
                --refPointLabel TSS \
                --yAxisLabel "Average signal" \
                --plotTitle "{wildcards.condition} {wildcards.assay} (averaged replicates)" \
                -out {output.pdf} \
                >> {log} 2>&1
            """
    
    # Gene body profiles - averaged (mean signal)
    rule profile_genes_mean:
        input:
            genes=os.path.join(OUT_DEEPTOOLS, "genes.bed"),
            bw=os.path.join("bigwig_averaged", "{condition}_{assay}_mean.bw")
        output:
            matrix=os.path.join("condition_profiles", "{condition}_{assay}_genes_mean.mat.gz"),
            pdf=os.path.join("condition_profiles", "{condition}_{assay}_genes_mean_profile.pdf")
        threads: THREADS
        log:
            os.path.join(LOG_DIR, "{condition}_{assay}_genes_mean_profile.log")
        shell:
            """
            mkdir -p condition_profiles
            
            computeMatrix scale-regions \
                -b 3000 -a 3000 \
                --regionBodyLength 5000 \
                --binSize 250 \
                -R {input.genes} \
                -S {input.bw} \
                --skipZeros --missingDataAsZero \
                -p {threads} \
                -o {output.matrix} \
                > {log} 2>&1
            
            plotProfile \
                -m {output.matrix} \
                --yAxisLabel "Average signal" \
                --plotTitle "{wildcards.condition} {wildcards.assay} gene-body (averaged replicates)" \
                -out {output.pdf} \
                >> {log} 2>&1
            """

