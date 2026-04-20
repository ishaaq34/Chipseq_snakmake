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
#   - For "trim" mode: trimmed FASTQ files in the configured trim directory
#   - For "bam" mode: BAMs in the configured pre-dedup alignment directory
# =====================================================

import os
import sys

# =====================================================
# Configuration and Report Setup
# =====================================================

# Global report configuration
report: "report/workflow.rst"

# Extract configuration
GENOME_FASTA = config["genome"]["fasta"]
GTF_FILE = config["genome"]["gtf"]
GENOME_INDEX = config["genome"]["index_prefix"]
EFFECTIVE_GENOME_SIZE = config["genome"]["effective_size"]

RAW_READS_DIR = config["directories"]["raw_reads"]
PREPROCESSING_DIR = config["directories"]["preprocessing"]
ALIGNMENT_DIR = config["directories"]["alignment"]
BAM_QC_DIR = config["directories"]["bam_qc"]
OUT_DEEPTOOLS = config["directories"]["deeptools"]
LOG_DIR = config["directories"]["logs"]
PEAK_ANALYSIS_DIR = config["directories"]["peak_analysis"]
REPORTS_DIR = config["directories"]["reports"]

TRIM_DIR = PREPROCESSING_DIR
OUT_ALIGN = os.path.join(ALIGNMENT_DIR, "1.1_pre_dedup_bam")
OUT_DEDUP = os.path.join(ALIGNMENT_DIR, "1.2_post_dedup_bam")
OUT_QC = os.path.join(BAM_QC_DIR, "2.1_general_bam_qc")
CHIP_QC_DIR = os.path.join(BAM_QC_DIR, "2.2_chipseq_specific_bam_qc")

SIGNAL_TRACKS_DIR = os.path.join(OUT_DEEPTOOLS, "3.1_signal_tracks")
PER_SAMPLE_BIGWIG_DIR = os.path.join(SIGNAL_TRACKS_DIR, "3.1.1_per_sample_bigwigs")
REFERENCE_REGIONS_DIR = os.path.join(SIGNAL_TRACKS_DIR, "3.1.2_reference_regions")
GROUPED_QC_DIR = os.path.join(OUT_DEEPTOOLS, "3.2_replicatewise_qc_profiles")
ENRICHMENT_ANALYSIS_DIR = os.path.join(OUT_DEEPTOOLS, "3.3_averaged_and_input_normalized_enrichment_outputs")
AVERAGED_BIGWIG_DIR = os.path.join(ENRICHMENT_ANALYSIS_DIR, "3.3.1_replicate_averaged_bigwigs")
NORMALIZED_BIGWIG_DIR = os.path.join(ENRICHMENT_ANALYSIS_DIR, "3.3.2_chip_over_input_log2_bigwigs")
ENRICHMENT_PROFILES_DIR = os.path.join(ENRICHMENT_ANALYSIS_DIR, "3.3.3_tss_gene_body_enrichment_profiles")

PEAK_CALL_DIR = os.path.join(PEAK_ANALYSIS_DIR, "4.1_peak_calling")
FRIP_DIR = os.path.join(PEAK_ANALYSIS_DIR, "4.2_frip_scores")
IDR_DIR = os.path.join(PEAK_ANALYSIS_DIR, "4.3_idr_reproducibility")
CONSENSUS_DIR = os.path.join(PEAK_ANALYSIS_DIR, "4.4_consensus_peaks")
HOMER_DIR = os.path.join(PEAK_ANALYSIS_DIR, "4.5_motif_enrichment")
ANALYSIS_SUMMARY_PATH = os.path.join(REPORTS_DIR, "analysis_summary.xlsx")

THREADS = config["threads"]
BIN_SIZE = config["deeptools"]["bin_size"]
SMOOTH_LENGTH = config["deeptools"]["smooth_length"]
MAPQ_THRESHOLD = config["samtools"]["mapq_threshold"]
SAM_FLAGS = config["samtools"]["sam_flags_exclude"]
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
    
    print(f"✓ Metadata validation passed: {len(df)} rows, {len(assays)} assays, {df['condition'].nunique()} conditions", file=sys.stderr)
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
    print(f"Using metadata TSV: {config['metadata_tsv']}", file=sys.stderr)
    metadata_df = validate_metadata_tsv(config["metadata_tsv"])
    MACS3_SAMPLES_CONFIG = derive_chip_input_pairs(metadata_df)
    print(f"✓ Auto-paired {len([s for s in MACS3_SAMPLES_CONFIG if MACS3_SAMPLES_CONFIG[s]['mark_type']])} ChIP samples with inputs", file=sys.stderr)
    
    # Extract sample IDs from metadata (no separate samples_file needed!)
    SAMPLES = list(metadata_df["raw_id"].unique())
    print(f"✓ Extracted {len(SAMPLES)} samples from metadata TSV", file=sys.stderr)
    
elif config.get("macs3_samples"):
    # Mode 1: Explicit (backward compatible)
    MACS3_SAMPLES_CONFIG = config["macs3_samples"]
    if MACS3_SAMPLES_CONFIG:
        print(f"✓ Using explicit sample configuration ({len(MACS3_SAMPLES_CONFIG)} samples)", file=sys.stderr)
    
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
# Get samples that need peak calling (exclude inputs)
MACS3_SAMPLES = [s for s in MACS3_SAMPLES_CONFIG.keys() 
                 if MACS3_SAMPLES_CONFIG[s].get("mark_type")]

# -----------------------------------------------------
# DeepTools QC Grouping Logic
# -----------------------------------------------------
def build_deeptools_groups():
    """
    Build sample groups for DeepTools QC based on config 'qc_grouping'.
    Returns dict: group_name -> [list of sample IDs]
    """
    modes = config.get("deeptools", {}).get("qc_grouping", ["all"])
    groups = {}

    # Mode 1: All samples
    if "all" in modes:
        groups["all_samples"] = SAMPLES

    # Mode 2: Per-condition
    if "condition" in modes and "metadata_tsv" in config:
        for cond in metadata_df["condition"].unique():
            samples = list(metadata_df[metadata_df["condition"] == cond]["raw_id"])
            if samples:
                groups[cond] = samples

    # Mode 3: Per-condition + Assay (always includes matched Input)
    if "condition_assay" in modes and "metadata_tsv" in config:
        for cond in metadata_df["condition"].unique():
            for assay in metadata_df[metadata_df["condition"] == cond]["assay"].unique():
                if assay.lower() == "input":
                    continue
                
                # ChIP samples for this condition+assay
                chip_samples = list(metadata_df[
                    (metadata_df["condition"] == cond) & 
                    (metadata_df["assay"] == assay)
                ]["raw_id"])
                
                # Matched Input samples for this condition
                input_samples = list(metadata_df[
                    (metadata_df["condition"] == cond) & 
                    (metadata_df["assay"].str.lower() == "input")
                ]["raw_id"])
                
                group_name = f"{cond}_{assay}"
                groups[group_name] = sorted(list(set(chip_samples + input_samples)))

    # Fallback if no metadata
    if not groups and SAMPLES:
        groups["all_samples"] = SAMPLES
        
    return groups

DEEPTOOLS_GROUPS = build_deeptools_groups()

def get_group_bams(wildcards):
    """Input function for DeepTools rules"""
    if wildcards.group not in DEEPTOOLS_GROUPS:
        return []
    return [os.path.join(OUT_DEDUP, f"{s}.dedup.bam") for s in DEEPTOOLS_GROUPS[wildcards.group]]

def get_group_labels(wildcards):
    """Labels for DeepTools plots"""
    if wildcards.group not in DEEPTOOLS_GROUPS:
        return ""
    return " ".join(DEEPTOOLS_GROUPS[wildcards.group])

def get_group_bigwigs(wildcards):
    """BigWig input function for grouped profile plots"""
    if wildcards.group not in DEEPTOOLS_GROUPS:
        return []
    return [os.path.join(PER_SAMPLE_BIGWIG_DIR, f"{s}.dedup.bw") for s in DEEPTOOLS_GROUPS[wildcards.group]]

def get_replicate_pairs_for_assay(condition, assay):
    """Get all pairwise replicate combinations for a condition and assay"""
    from itertools import combinations
    
    samples = [s for s, info in MACS3_SAMPLES_CONFIG.items() 
               if info.get("mark_type") and info.get("condition") == condition and info.get("assay") == assay]
    
    reps = sorted([MACS3_SAMPLES_CONFIG[s].get("replicate") for s in samples])
    return list(combinations(reps, 2))

def get_sample_for_replicate_assay(condition, assay, replicate):
    """Get sample name for a given condition, assay and replicate number"""
    for sample, info in MACS3_SAMPLES_CONFIG.items():
        if (info.get("condition") == condition and 
            info.get("assay") == assay and 
            str(info.get("replicate")) == str(replicate)):
            if info.get("mark_type"):  # ChIP sample (not input)
                return sample
    return None

# Helper functions for peak type analysis
def get_peak_types_for_sample(sample):
    """
    Determine which peak types to call for a sample.
    Returns: list of peak types ["narrow"], ["broad"], or ["narrow", "broad"]
    """
    mark_type = MACS3_SAMPLES_CONFIG.get(sample, {}).get("mark_type", "narrow")
    
    if mark_type == "both":
        return ["narrow", "broad"]
    elif mark_type == "broad":
        return ["broad"]
    else:  # narrow or empty
        return ["narrow"]

def get_all_peak_types_for_condition_assay(condition, assay):
    """
    Get all peak types for a given condition/assay combination.
    Combines peak types from all replicates.
    """
    peak_types = set()
    for sample, info in MACS3_SAMPLES_CONFIG.items():
        if info.get("condition") == condition and info.get("assay") == assay:
            peak_types.update(get_peak_types_for_sample(sample))
    return sorted(list(peak_types))

def get_samples_for_condition_assay(condition, assay):
    """
    Get all sample names for a given condition/assay combination.
    """
    samples = []
    for sample, info in MACS3_SAMPLES_CONFIG.items():
        if info.get("condition") == condition and info.get("assay") == assay:
            samples.append(sample)
    return samples


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
os.makedirs(CHIP_QC_DIR, exist_ok=True)
os.makedirs(OUT_DEEPTOOLS, exist_ok=True)
os.makedirs(LOG_DIR, exist_ok=True)
os.makedirs(PEAK_ANALYSIS_DIR, exist_ok=True)
os.makedirs(REPORTS_DIR, exist_ok=True)
os.makedirs(SIGNAL_TRACKS_DIR, exist_ok=True)
os.makedirs(ENRICHMENT_ANALYSIS_DIR, exist_ok=True)
os.makedirs(PER_SAMPLE_BIGWIG_DIR, exist_ok=True)
os.makedirs(REFERENCE_REGIONS_DIR, exist_ok=True)
os.makedirs(GROUPED_QC_DIR, exist_ok=True)
os.makedirs(AVERAGED_BIGWIG_DIR, exist_ok=True)
os.makedirs(NORMALIZED_BIGWIG_DIR, exist_ok=True)
os.makedirs(ENRICHMENT_PROFILES_DIR, exist_ok=True)
os.makedirs(PEAK_CALL_DIR, exist_ok=True)
os.makedirs(FRIP_DIR, exist_ok=True)
os.makedirs(IDR_DIR, exist_ok=True)
os.makedirs(CONSENSUS_DIR, exist_ok=True)
os.makedirs(HOMER_DIR, exist_ok=True)

# =====================================================
# Target Outputs
# =====================================================
ALL_BIGWIGS = expand(os.path.join(PER_SAMPLE_BIGWIG_DIR, "{sample}.dedup.bw"), sample=SAMPLES)
ALL_FLAGS = expand(os.path.join(OUT_QC, "{sample}.flagstat.txt"), sample=SAMPLES)
ALL_PBC = expand(os.path.join(CHIP_QC_DIR, "{sample}.pbc.txt"), sample=SAMPLES)
ALL_SPP = expand([
    os.path.join(CHIP_QC_DIR, "{sample}_spp.qc.txt"),
    os.path.join(CHIP_QC_DIR, "{sample}_avp.pdf"),
], sample=SAMPLES) if RUN_SPP else []

# MultiQC reports
MULTIQC_ALIGNMENT = os.path.join(OUT_QC, "alignment_qc_report.html")

ALL_DEEPTOOLS_QC = []
for grp in DEEPTOOLS_GROUPS:
    ALL_DEEPTOOLS_QC.extend([
        os.path.join(GROUPED_QC_DIR, grp, "all_samples_bins.npz"),
        os.path.join(GROUPED_QC_DIR, grp, "correlation_heatmap.pdf"),
        os.path.join(GROUPED_QC_DIR, grp, "correlation_scatterplot.pdf"),
        os.path.join(GROUPED_QC_DIR, grp, "pca.pdf"),
        os.path.join(GROUPED_QC_DIR, grp, "fingerprints.pdf"),
        os.path.join(GROUPED_QC_DIR, grp, "coverage_histogram.pdf"),
    ])
ALL_DEEPTOOLS_QC = [x for x in ALL_DEEPTOOLS_QC if x]  # Remove None values

ALL_PROFILE_OUTPUTS = [
    os.path.join(GROUPED_QC_DIR, "matrix_all_bw_TSS.gz"),
    os.path.join(GROUPED_QC_DIR, "all_samples_TSS_profile.pdf"),
    os.path.join(GROUPED_QC_DIR, "matrix_all_bw_scalar.gz"),
    os.path.join(GROUPED_QC_DIR, "all_samples_scalar_profile.pdf")
]

# Grouped profile plots (TSS + gene body per qc_grouping)
for grp in DEEPTOOLS_GROUPS:
    ALL_PROFILE_OUTPUTS.extend([
        os.path.join(GROUPED_QC_DIR, grp, "TSS_profile.pdf"),
        os.path.join(GROUPED_QC_DIR, grp, "gene_body_profile.pdf"),
    ])

# MACS3 peak calling outputs (conditional on configuration)
# MACS3 peak calling outputs (hierarchical structure)
ALL_MACS3_PEAKS = []
if MACS3_SAMPLES:
    for sample in MACS3_SAMPLES:
        sample_info = MACS3_SAMPLES_CONFIG[sample]
        condition = sample_info.get("condition")
        assay = sample_info.get("assay")
        
        # Get peak types for this sample (narrow, broad, or both)
        peak_types = get_peak_types_for_sample(sample)
        
        for peak_type in peak_types:
            if peak_type == "narrow":
                ALL_MACS3_PEAKS.extend([
                    os.path.join(PEAK_CALL_DIR, condition, assay, "narrow", f"{sample}_peaks.narrowPeak"),
                    os.path.join(PEAK_CALL_DIR, condition, assay, "narrow", f"{sample}_summits.bed")
                ])
            else:  # broad
                ALL_MACS3_PEAKS.extend([
                    os.path.join(PEAK_CALL_DIR, condition, assay, "broad", f"{sample}_peaks.broadPeak"),
                    os.path.join(PEAK_CALL_DIR, condition, assay, "broad", f"{sample}_peaks.gappedPeak")
                ])


# FRiP calculation outputs (hierarchical structure)
ALL_FRIP = []
if MACS3_SAMPLES:
    for sample in MACS3_SAMPLES:
        sample_info = MACS3_SAMPLES_CONFIG[sample]
        condition = sample_info.get("condition")
        assay = sample_info.get("assay")
        
        # Get peak types for this sample
        peak_types = get_peak_types_for_sample(sample)
        
        for peak_type in peak_types:
            ALL_FRIP.append(os.path.join(FRIP_DIR, condition, assay, peak_type, f"{sample}.frip.txt"))
    # Aggregate table - skip for now, will update later
    # ALL_FRIP.append(os.path.join(FRIP_DIR, "all_samples_frip.tsv"))
# DISABLED - requires datavzrd config:     ALL_FRIP.append(os.path.join(FRIP_DIR, "interactive_table", "index.html"))

# IDR analysis outputs (conditional on metadata TSV with ≥2 replicates)
ALL_IDR = []
ALL_CONSENSUS_PEAKS = []
if MACS3_SAMPLES_CONFIG and config.get("metadata_tsv"):
    # Group samples by condition
    from itertools import combinations
    
    # Group ChIP samples by (condition, assay) for IDR and Consensus
    condition_assay_chip = {}
    for sample, info in MACS3_SAMPLES_CONFIG.items():
        if info.get("mark_type"):  # Only ChIP samples
            cond = info.get("condition", "default")
            assay = info.get("assay", "default")
            key = (cond, assay)
            if key not in condition_assay_chip:
                condition_assay_chip[key] = []
            condition_assay_chip[key].append(sample)
    
    # Generate pairwise IDR comparisons and consensus peaks per peak type
    for (condition, assay), samples in condition_assay_chip.items():

        peak_types = get_all_peak_types_for_condition_assay(condition, assay)
        for peak_type in peak_types:
            if len(samples) >= 2:
                for sample1, sample2 in combinations(sorted(samples), 2):
                    rep1 = MACS3_SAMPLES_CONFIG[sample1]["replicate"]
                    rep2 = MACS3_SAMPLES_CONFIG[sample2]["replicate"]
                    
                    ALL_IDR.append(os.path.join(IDR_DIR, condition, assay, peak_type, f"rep{rep1}_vs_rep{rep2}_idr.txt"))
                    ALL_IDR.append(os.path.join(IDR_DIR, condition, assay, peak_type, f"rep{rep1}_vs_rep{rep2}_idr.txt.png"))
            
            # Consensus peaks per peak type
            ALL_CONSENSUS_PEAKS.append(os.path.join(CONSENSUS_DIR, condition, assay, peak_type, "consensus.bed"))
            ALL_CONSENSUS_PEAKS.append(os.path.join(CONSENSUS_DIR, condition, assay, peak_type, "consensus_frip.txt"))


# HOMER motif analysis outputs (conditional on run_homer flag)
ALL_HOMER = []
if RUN_HOMER and ALL_CONSENSUS_PEAKS:
    for (condition, assay), samples in condition_assay_chip.items():
        peak_types = get_all_peak_types_for_condition_assay(condition, assay)
        for peak_type in peak_types:
            ALL_HOMER.append(os.path.join(HOMER_DIR, condition, assay, peak_type, "homerResults.html"))



# Averaged BigWigs per condition/assay (conditional on metadata TSV)
ALL_AVERAGED_BIGWIGS = []
if MACS3_SAMPLES_CONFIG and config.get("metadata_tsv"):
    # Group ALL samples by (condition, assay) for averaging
    condition_assay_all = {}
    for sample, info in MACS3_SAMPLES_CONFIG.items():
        cond = info.get("condition", "default")
        assay = info.get("assay", "unknown")
        key = (cond, assay)
        
        if key not in condition_assay_all:
            condition_assay_all[key] = []
        condition_assay_all[key].append(sample)
    
    # Only create averaged BigWigs for groups with ≥2 replicates
    for (condition, assay), samples in condition_assay_all.items():
        if len(samples) >= 2:
            ALL_AVERAGED_BIGWIGS.append(os.path.join(AVERAGED_BIGWIG_DIR, condition, f"{assay}_mean.bw"))



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
            ALL_NORMALIZED_BIGWIGS.append(os.path.join(NORMALIZED_BIGWIG_DIR, condition, f"{assay}_log2_over_input.bw"))


# Profile plots for averaged/normalized BigWigs per condition
ALL_CONDITION_PROFILES = []
if MACS3_SAMPLES_CONFIG and config.get("metadata_tsv"):
    for condition, assays in chip_conditions.items():
        for assay in assays:
            # Normalized BigWig profiles (log2 ratio)
            ALL_CONDITION_PROFILES.append(os.path.join(ENRICHMENT_PROFILES_DIR, condition, assay, "TSS_log2_profile.pdf"))
            ALL_CONDITION_PROFILES.append(os.path.join(ENRICHMENT_PROFILES_DIR, condition, assay, "genes_log2_profile.pdf"))
            # Averaged BigWig profiles (mean signal)
            ALL_CONDITION_PROFILES.append(os.path.join(ENRICHMENT_PROFILES_DIR, condition, assay, "TSS_mean_profile.pdf"))
            ALL_CONDITION_PROFILES.append(os.path.join(ENRICHMENT_PROFILES_DIR, condition, assay, "genes_mean_profile.pdf"))


# =====================================================
# Main Rule
# =====================================================
rule all:
    input:
        ALL_BIGWIGS,
        ALL_FLAGS,
        ALL_PBC,
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
        ALL_CONDITION_PROFILES,
        os.path.join(CHIP_QC_DIR, "chipseq_qc_summary.tsv"),  # Combined ChIP QC summary
        ANALYSIS_SUMMARY_PATH,  # Comprehensive Excel summary
        "docs/results/latest_run.md"  # Automated MkDocs results summary

# =====================================================
# 1. Quality Filtering and Adapter Trimming (Trim_galore)
# =====================================================
if START_FROM in ["raw", "fastq"]:
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
if START_FROM in ["raw", "fastq", "trim"]:
    rule build_index:
        input:
            GENOME_FASTA
        output:
            expand(
                GENOME_INDEX + ".{ext}.bt2",
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
# Use conditional rule definition based on sequencing type (PE or SE)
# Only ONE rule will be defined at parse-time based on config

if START_FROM in ["raw", "fastq", "trim"]:
    if SEQ_TYPE == "paired":
        rule align_sort:
            input:
                r1=os.path.join(TRIM_DIR, "{sample}_R1_val_1.fq.gz"),
                r2=os.path.join(TRIM_DIR, "{sample}_R2_val_2.fq.gz"),
                index=expand(
                    os.path.join(os.path.dirname(GENOME_INDEX) or ".", os.path.basename(GENOME_INDEX) + ".{ext}.bt2"),
                    ext=["1", "2", "3", "4", "rev.1", "rev.2"]
                )
            output:
                bam=os.path.join(OUT_ALIGN, "{sample}.sorted.bam"),
                bai=os.path.join(OUT_ALIGN, "{sample}.sorted.bam.bai")
            params:
                index=GENOME_INDEX,
                sensitivity=config["bowtie2"].get("sensitivity", "--sensitive"),
                no_unal="--no-unal" if config["bowtie2"].get("no_unal", True) else "",
                no_mixed="--no-mixed" if config["bowtie2"].get("no_mixed", True) else "",
                no_discordant="--no-discordant" if config["bowtie2"].get("no_discordant", True) else "",
                mapq=config["samtools"].get("mapq_threshold", 30),
                flags=config["samtools"].get("sam_flags_exclude", 768)
            threads: THREADS
            resources:
                mem_mb=config["resources"]["align_mem_mb"]
            log:
                bowtie2=os.path.join(LOG_DIR, "{sample}_bowtie2.log"),
                samtools_sort=os.path.join(LOG_DIR, "{sample}_samtools_sort.log")
            shell:
                """
                bowtie2 \\
                  -p {threads} \\
                  -x {params.index} \\
                  -1 {input.r1} \\
                  -2 {input.r2} \\
                  {params.sensitivity} \\
                  {params.no_mixed} \\
                  {params.no_discordant} \\
                  {params.no_unal} \\
                  2> {log.bowtie2} \\
                | samtools view -@ {threads} -bS -q {params.mapq} -F {params.flags} - \\
                | samtools sort -@ {threads} -o {output.bam} - 2> {log.samtools_sort}
                
                samtools index -@ {threads} {output.bam}
                """
    
    else:  # Single-end
        rule align_sort:
            input:
                r1=os.path.join(TRIM_DIR, "{sample}_trimmed.fq.gz"),
                index=expand(
                    os.path.join(os.path.dirname(GENOME_INDEX) or ".", os.path.basename(GENOME_INDEX) + ".{ext}.bt2"),
                    ext=["1", "2", "3", "4", "rev.1", "rev.2"]
                )
            output:
                bam=os.path.join(OUT_ALIGN, "{sample}.sorted.bam"),
                bai=os.path.join(OUT_ALIGN, "{sample}.sorted.bam.bai")
            params:
                index=GENOME_INDEX,
                sensitivity=config["bowtie2"].get("sensitivity", "--very-sensitive"),
                no_unal="--no-unal" if config["bowtie2"].get("no_unal", True) else "",
                mapq=config["samtools"].get("mapq_threshold", 30),
                flags=config["samtools"].get("sam_flags_exclude", 768)
            threads: THREADS
            resources:
                mem_mb=config["resources"]["align_mem_mb"]
            log:
                bowtie2=os.path.join(LOG_DIR, "{sample}_bowtie2.log"),
                samtools_sort=os.path.join(LOG_DIR, "{sample}_samtools_sort.log")
            shell:
                """
                bowtie2 \\
                  -p {threads} \\
                  -x {params.index} \\
                  -U {input.r1} \\
                  {params.sensitivity} \\
                  {params.no_unal} \\
                  2> {log.bowtie2} \\
                | samtools view -@ {threads} -bS -q {params.mapq} -F {params.flags} - \\
                | samtools sort -@ {threads} -o {output.bam} - 2> {log.samtools_sort}
                
                samtools index -@ {threads} {output.bam}
                """

if START_FROM == "bam":
    rule sort_input_bam:
        input:
            bam=os.path.join(OUT_ALIGN, "{sample}.bam")
        output:
            bam=os.path.join(OUT_ALIGN, "{sample}.sorted.bam"),
            bai=os.path.join(OUT_ALIGN, "{sample}.sorted.bam.bai")
        threads: THREADS
        resources:
            mem_mb=config["resources"]["align_mem_mb"]
        log:
            os.path.join(LOG_DIR, "{sample}_sort_input_bam.log")
        shell:
            """
            samtools sort -@ {threads} -o {output.bam} {input.bam} > {log} 2>&1
            samtools index -@ {threads} {output.bam} >> {log} 2>&1
            """
# Interactive FRiP table with Datavzrd
rule view_frip_with_datavzrd:
    input:
        config="datavzrd/frip_summary.yaml",
        table=os.path.join(FRIP_DIR, "all_samples_frip.tsv")
    output:
        frip_report=report(
            directory(os.path.join(FRIP_DIR, "interactive_table")),
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
    """
    Perform duplicate marking and optional removal using samtools markdup.
    The behavior is controlled by config["dedup"]["remove_duplicates"]:
    - if true: duplicates are removed (-r flag is added)
    - if false: duplicates are only marked (no -r flag)
    """
    input:
        bam=os.path.join(OUT_ALIGN, "{sample}.sorted.bam")
    output:
        bam=os.path.join(OUT_DEDUP, "{sample}.dedup.bam")
    resources:
        mem_mb=config["resources"].get("markdup_mem_mb", 4000)
    params:
        # Robustly extract the remove_duplicates flag
        remove_dup_flag=lambda wildcards: "-r" if config.get("dedup", {}).get("remove_duplicates", True) else ""
    log:
        os.path.join(LOG_DIR, "{sample}_markdup.log")
    threads: THREADS
    shell:
        """
        mkdir -p $(dirname {output.bam})
        
        # Determine if we should remove or just mark duplicates based on config
        FLAG="{params.remove_dup_flag}"
        
        if [ "{SEQ_TYPE}" == "paired" ]; then
            # Paired-end deduplication using piped workflow (Strict HTSlib adherence)
            # Reference: https://www.htslib.org/algorithms/duplicate.html
            # Workflow: Collate (name group) -> Fixmate (add mate scores) -> Sort (position) -> Markdup
            samtools collate -@ {threads} -O -u {input.bam} 2> {log} | \
            samtools fixmate -@ {threads} -m -u - - 2>> {log} | \
            samtools sort -@ {threads} -u - 2>> {log} | \
            samtools markdup -@ {threads} $FLAG - {output.bam} 2>> {log}
        else
            # Single-end deduplication
            # Reference: https://www.htslib.org/algorithms/duplicate.html
            # Requirement: Input must be coordinate sorted (satisfied by upstream 'align_sort' rule)
            samtools markdup -@ {threads} $FLAG {input.bam} {output.bam} 2>> {log}
        fi
        
        samtools index -@ {threads} {output.bam}
        """

# =====================================================
# 5. Quality Control Metrics
# =====================================================
rule qc_flagstat_stats:
    input:
        bam=os.path.join(OUT_DEDUP, "{sample}.dedup.bam")
    output:
        flagstat=os.path.join(OUT_QC, "{sample}.flagstat.txt"),
        stats=os.path.join(OUT_QC, "{sample}.stats.txt")
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
        """
# =====================================================
# 5b. ChIP-seq Library Complexity (PBC/NRF)
# =====================================================
rule qc_pbc:
    input:
        bam=os.path.join(OUT_DEDUP, "{sample}.dedup.bam")
    output:
        read5=temp(os.path.join(CHIP_QC_DIR, "{sample}.read5.bed")),
        pbc=os.path.join(CHIP_QC_DIR, "{sample}.pbc.txt")
    threads: 1
    resources:
        mem_mb=config["resources"]["default_mem_mb"]
    log:
        os.path.join(LOG_DIR, "{sample}_pbc.log")
    shell:
        """
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
        alignment_qc_report=report(
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
            avp=os.path.join(CHIP_QC_DIR, "{sample}_avp.pdf"),
            sppqc=os.path.join(CHIP_QC_DIR, "{sample}_spp.qc.txt")
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
# 7b. Combined ChIP-seq QC Summary (PBC + SPP)
# =====================================================
rule chipseq_qc_summary:
    input:
        pbc=expand(os.path.join(CHIP_QC_DIR, "{sample}.pbc.txt"), sample=SAMPLES),
        spp=expand(os.path.join(CHIP_QC_DIR, "{sample}_spp.qc.txt"), sample=SAMPLES) if RUN_SPP else [],
        metadata=config.get("metadata_tsv", "metadata_from_new_names.tsv")
    output:
        tsv=os.path.join(CHIP_QC_DIR, "chipseq_qc_summary.tsv"),
        nsc_png=os.path.join(CHIP_QC_DIR, "NSC_plot.png"),
        rsc_png=os.path.join(CHIP_QC_DIR, "RSC_plot.png"),
        nrf_png=os.path.join(CHIP_QC_DIR, "NRF_plot.png"),
        pbc1_png=os.path.join(CHIP_QC_DIR, "PBC1_plot.png"),
        pbc2_png=os.path.join(CHIP_QC_DIR, "PBC2_plot.png"),
    log:
        os.path.join(LOG_DIR, "chipseq_qc_summary.log")
    run:
        import pandas as pd
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        import numpy as np
        import os, glob, re

        # ── 1. Read metadata ──
        meta = pd.read_csv(input.metadata, sep="\t")
        meta = meta.rename(columns={"raw_id": "Sample"})

        # ── 2. Parse PBC files ──
        pbc_rows = []
        for f in input.pbc:
            sample = os.path.basename(f).replace(".pbc.txt", "")
            with open(f) as fh:
                line = fh.read().strip()
            vals = {}
            for part in line.split("\t"):
                k, v = part.split("=")
                vals[k] = float(v) if v != "Inf" else float("inf")
            vals["Sample"] = sample
            pbc_rows.append(vals)
        pbc_df = pd.DataFrame(pbc_rows)

        # ── 3. Parse SPP files ──
        spp_rows = []
        spp_cols = ["Filename", "NumReads", "EstFragLen", "Corr_EstFragLen",
                     "PhantomPeak", "Corr_PhantomPeak", "ArgMin_Corr",
                     "Min_Corr", "NSC", "RSC", "QualityTag"]
        for f in sorted(glob.glob(os.path.join(CHIP_QC_DIR, "*_spp.qc.txt"))):
            sample = os.path.basename(f).replace("_spp.qc.txt", "")
            with open(f) as fh:
                line = fh.read().strip()
            parts = line.split("\t")
            row = {"Sample": sample}
            for i, col in enumerate(spp_cols):
                if i < len(parts):
                    if col in ["NumReads", "PhantomPeak", "ArgMin_Corr", "QualityTag"]:
                        try:
                            row[col] = int(parts[i])
                        except:
                            row[col] = parts[i]
                    elif col in ["NSC", "RSC", "Corr_PhantomPeak", "Min_Corr"]:
                        try:
                            row[col] = float(parts[i])
                        except:
                            row[col] = parts[i]
                    elif col == "EstFragLen":
                        # Take first value if comma-separated
                        row[col] = int(parts[i].split(",")[0])
                        row["EstFragLen_All"] = parts[i]
                    elif col == "Corr_EstFragLen":
                        row[col] = float(parts[i].split(",")[0])
                    else:
                        row[col] = parts[i]
            spp_rows.append(row)
        spp_df = pd.DataFrame(spp_rows)
        if "Filename" in spp_df.columns:
            spp_df = spp_df.drop(columns=["Filename"])

        # ── 4. Merge all ──
        df = meta.merge(pbc_df, on="Sample", how="left")
        if len(spp_df) > 0:
            df = df.merge(spp_df, on="Sample", how="left")

        # ── 5. ENCODE QC flags ──
        if "NSC" in df.columns:
            df["NSC_Pass"] = df["NSC"].apply(lambda x: "PASS" if pd.notna(x) and x > 1.05 else "FAIL")
        if "RSC" in df.columns:
            df["RSC_Pass"] = df["RSC"].apply(lambda x: "PASS" if pd.notna(x) and x > 0.8 else "FAIL")
        df["NRF_Pass"] = df["NRF"].apply(lambda x: "PASS" if pd.notna(x) and x > 0.8 else "FAIL")

        # Sort for display — "input" always last, others alphabetical
        unique_assays = sorted([a for a in df["assay"].unique() if a.lower() != "input"])
        if "input" in df["assay"].str.lower().values:
            input_name = df.loc[df["assay"].str.lower() == "input", "assay"].iloc[0]
            unique_assays.append(input_name)
        assay_order = {a: i for i, a in enumerate(unique_assays)}
        df["_sort"] = df["assay"].map(assay_order).fillna(99)
        df = df.sort_values(["condition", "_sort", "replicate"]).drop(columns=["_sort"])

        # ── 6. Save summary TSV ──
        out_cols = ["Sample", "condition", "assay", "replicate"]
        if "NumReads" in df.columns:
            out_cols += ["NumReads", "EstFragLen", "NSC", "RSC", "QualityTag"]
        out_cols += ["NRF", "PBC1", "PBC2"]
        if "NSC_Pass" in df.columns:
            out_cols += ["NSC_Pass", "RSC_Pass"]
        out_cols += ["NRF_Pass"]
        existing = [c for c in out_cols if c in df.columns]
        df[existing].to_csv(output.tsv, sep="\t", index=False)
        print(f"Summary saved: {output.tsv}")

        # ── 7. Plots ──
        # Auto-generate color palette from metadata assays
        # Palette is configurable via config["plotting"]["palette"] (default: Dark2)
        # Any matplotlib/RColorBrewer palette works: Dark2, Set1, tab10, viridis, Paired, etc.
        palette_name = config.get("plotting", {}).get("palette", "Dark2")
        all_assays = unique_assays  # already sorted with input last
        cmap = plt.cm.get_cmap(palette_name, max(len(all_assays), 3))
        assay_colors = {}
        for i, assay in enumerate(all_assays):
            if assay.lower() == "input":
                assay_colors[assay] = "#999999"  # neutral gray for input controls
            else:
                assay_colors[assay] = matplotlib.colors.rgb2hex(cmap(i))

        # Auto-generate markers from metadata conditions
        unique_conditions = sorted(df["condition"].unique())
        marker_pool = ["o", "s", "D", "^" , "v", "P", "X", "*"]
        cond_markers = {cond: marker_pool[i % len(marker_pool)] for i, cond in enumerate(unique_conditions)}

        # Helper: create grouped bar plot and save as individual figure
        def make_bar_plot(df, metric, threshold, ylabel, title, out_path):
            fig, ax = plt.subplots(figsize=(10, 6))
            # Group by condition + assay
            groups = df.groupby(["condition", "assay"], sort=False)
            labels = []
            positions = []
            colors_list = []
            values = []
            pos = 0
            prev_cond = None
            for (cond, assay), grp in groups:
                if prev_cond is not None and cond != prev_cond:
                    pos += 0.5  # gap between conditions
                for _, row in grp.iterrows():
                    labels.append(f"{row['Sample']}")
                    positions.append(pos)
                    colors_list.append(assay_colors.get(assay, "#999999"))
                    val = row[metric] if pd.notna(row[metric]) else 0
                    # Cap Inf values for display
                    if val == float("inf"):
                        val = 0
                    values.append(val)
                    pos += 1
                pos += 0.3  # gap between assays
                prev_cond = cond

            bars = ax.bar(positions, values, color=colors_list, width=0.8, edgecolor="white", linewidth=0.5)
            ax.axhline(y=threshold, color="#2D2D2D", linestyle="--", linewidth=1.5, alpha=0.7,
                       label=f"ENCODE threshold ({threshold})")
            ax.set_xticks(positions)
            ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=9)
            ax.set_ylabel(ylabel, fontsize=12)
            ax.set_title(title, fontsize=14, fontweight="bold")
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

            # Color bars that fail
            for bar, val in zip(bars, values):
                if val < threshold:
                    bar.set_alpha(0.5)
                    bar.set_hatch("//")

            # Build legend with assay colors + condition markers + threshold
            legend_handles = [plt.Line2D([0], [0], color="#2D2D2D", linestyle="--",
                              linewidth=1.5, label=f"ENCODE threshold ({threshold})")]
            for assay, color in assay_colors.items():
                legend_handles.append(mpatches.Patch(color=color, label=assay))
            ax.legend(handles=legend_handles, fontsize=9, loc="upper right")

            plt.tight_layout()
            plt.savefig(out_path, dpi=600, bbox_inches="tight")
            plt.close()
            print(f"Saved: {out_path}")

        # Generate individual plots
        if "NSC" in df.columns:
            make_bar_plot(df, "NSC", 1.05,
                         "Normalized Strand Coefficient",
                         "NSC (Signal-to-Noise)", output.nsc_png)

        if "RSC" in df.columns:
            make_bar_plot(df, "RSC", 0.8,
                         "Relative Strand Coefficient",
                         "RSC (Signal Enrichment)", output.rsc_png)

        make_bar_plot(df, "NRF", 0.8,
                     "Non-Redundant Fraction",
                     "NRF (Library Complexity)", output.nrf_png)

        make_bar_plot(df, "PBC1", 0.8,
                     "PCR Bottleneck Coefficient 1",
                     "PBC1 (PCR Bottleneck)", output.pbc1_png)

        make_bar_plot(df, "PBC2", 1.0,
                     "PCR Bottleneck Coefficient 2",
                     "PBC2 (PCR Bottleneck)", output.pbc2_png)

        print("All individual QC plots saved.")

# =====================================================
# 8. DeepTools QC (Multi-sample)
# =====================================================
rule grouped_multiBamSummary:
    input:
        bams=get_group_bams
    output:
        npz=os.path.join(GROUPED_QC_DIR, "{group}", "all_samples_bins.npz"),
        raw=os.path.join(GROUPED_QC_DIR, "{group}", "all_samples_bins.tab")
    params:
        labels=get_group_labels
    threads: THREADS
    resources:
        mem_mb=config["resources"]["deeptools_mem_mb"]
    log:
        os.path.join(LOG_DIR, "{group}_multiBamSummary.log")
    shell:
        """
        multiBamSummary bins \\
            -b {input.bams} \\
            --labels {params.labels} \\
            -o {output.npz} \\
            --outRawCounts {output.raw} \\
            --numberOfProcessors {threads} \\
            > {log} 2>&1
        """

rule grouped_plotCorrelation:
    input:
        npz=os.path.join(GROUPED_QC_DIR, "{group}", "all_samples_bins.npz")
    output:
        heatmap=report(
            os.path.join(GROUPED_QC_DIR, "{group}", "correlation_heatmap.pdf"),
            caption="report/correlation_plot.rst",
            category="QC Plots",
            labels={"plot": "correlation_heatmap", "group": "{group}"}
        ),
        scatter=report(
            os.path.join(GROUPED_QC_DIR, "{group}", "correlation_scatterplot.pdf"),
            caption="report/correlation_plot.rst",
            category="QC Plots",
            labels={"plot": "correlation_scatter", "group": "{group}"}
        ),
        matrix=os.path.join(GROUPED_QC_DIR, "{group}", "correlation_matrix.tab")
    params:
        dpi=config.get("deeptools", {}).get("plots", {}).get("dpi", 200),
        width=config.get("deeptools", {}).get("plots", {}).get("correlation", {}).get("plot_width", 11),
        height=config.get("deeptools", {}).get("plots", {}).get("correlation", {}).get("plot_height", 9.5)
    log:
        os.path.join(LOG_DIR, "{group}_plotCorrelation.log")
    shell:
        """
        # Generate heatmap
        plotCorrelation \\
            -in {input.npz} \\
            --corMethod spearman \\
            --skipZeros \\
            --whatToPlot heatmap \\
            --colorMap RdYlBu \\
            --plotNumbers \\
            --plotTitle "Spearman Correlation Heatmap ({wildcards.group})" \\
            --plotWidth {params.width} \\
            --plotHeight {params.height} \\
            -o {output.heatmap} \\
            --outFileCorMatrix {output.matrix} \\
            > {log} 2>&1
        
        # Generate scatterplot
        plotCorrelation \\
            -in {input.npz} \\
            --corMethod spearman \\
            --skipZeros \\
            --whatToPlot scatterplot \\
            --plotTitle "Spearman Correlation Scatterplot ({wildcards.group})" \\
            -o {output.scatter} \\
            >> {log} 2>&1
        """
rule grouped_plotPCA:
    input:
        npz=os.path.join(GROUPED_QC_DIR, "{group}", "all_samples_bins.npz")
    output:
        pca_plot=os.path.join(GROUPED_QC_DIR, "{group}", "pca.pdf"),
        tab=os.path.join(GROUPED_QC_DIR, "{group}", "PCA.tab")
    params:
        dpi=config.get("deeptools", {}).get("plots", {}).get("dpi", 200),
        width=config.get("deeptools", {}).get("plots", {}).get("pca", {}).get("plot_width", 14),
        height=config.get("deeptools", {}).get("plots", {}).get("pca", {}).get("plot_height", 12)
    log:
        os.path.join(LOG_DIR, "{group}_plotPCA.log")
    shell:
        """
        plotPCA \\
            -in {input.npz} \\
            -o {output.pca_plot} \\
            -T "PCA of Binned Coverage ({wildcards.group})" \\
            --transpose \\
            --plotWidth {params.width} \\
            --plotHeight {params.height} \\
            --outFileNameData {output.tab} \\
            > {log} 2>&1
        """
rule grouped_plotFingerprint:
    input:
        bams=get_group_bams
    output:
        fingerprints=report(
            os.path.join(GROUPED_QC_DIR, "{group}", "fingerprints.pdf"),
            caption="report/pca_plot.rst",
            category="QC Plots",
            labels={"plot": "fingerprint", "group": "{group}"}
        ),
        counts=os.path.join(GROUPED_QC_DIR, "{group}", "fingerprints.tab")
    params:
        labels=get_group_labels
    threads: THREADS
    resources:
        mem_mb=config["resources"]["deeptools_mem_mb"]
    log:
        os.path.join(LOG_DIR, "{group}_plotFingerprint.log")
    shell:
        """
        plotFingerprint \\
            -b {input.bams} \\
            --labels {params.labels} \\
            --skipZeros \\
            --numberOfSamples 50000 \\
            -T "Fingerprints ({wildcards.group})" \\
            --plotFile {output.fingerprints} \\
            --outRawCounts {output.counts} \\
            --numberOfProcessors {threads} \\
            > {log} 2>&1
        """

rule grouped_plotCoverage:
    input:
        bams=get_group_bams
    output:
        plot=os.path.join(GROUPED_QC_DIR, "{group}", "coverage_histogram.pdf"),
        counts=os.path.join(GROUPED_QC_DIR, "{group}", "coverage_counts.txt")
    params:
        labels=get_group_labels
    threads: THREADS
    log:
        os.path.join(LOG_DIR, "{group}_plotCoverage.log")
    shell:
        """
        plotCoverage \\
            -b {input.bams} \\
            --plotFile {output.plot} \\
            --labels {params.labels} \\
            --outRawCounts {output.counts} \\
            --numberOfProcessors {threads} \\
            > {log} 2>&1
        """
# =====================================================
# 10. BigWig Coverage Tracks
# =====================================================
rule bamCoverage:
    input:
        bam=os.path.join(OUT_DEDUP, "{sample}.dedup.bam")
    output:
        bw=os.path.join(PER_SAMPLE_BIGWIG_DIR, "{sample}.dedup.bw")
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
        tss=os.path.join(REFERENCE_REGIONS_DIR, "tss.bed"),
        genes=os.path.join(REFERENCE_REGIONS_DIR, "genes.bed")
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
        tss=os.path.join(REFERENCE_REGIONS_DIR, "tss.bed"),
        bws=expand(os.path.join(PER_SAMPLE_BIGWIG_DIR, "{sample}.dedup.bw"), sample=SAMPLES)
    output:
        matrix=os.path.join(GROUPED_QC_DIR, "matrix_all_bw_TSS.gz")
    params:
        upstream=config["deeptools"]["tss_upstream"],
        downstream=config["deeptools"]["tss_downstream"]
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
        matrix=os.path.join(GROUPED_QC_DIR, "matrix_all_bw_TSS.gz")
    output:
        profile=os.path.join(GROUPED_QC_DIR, "all_samples_TSS_profile.pdf"),
        data=os.path.join(GROUPED_QC_DIR, "all_samples_TSS_profile.tab")
    log:
        os.path.join(LOG_DIR, "plotProfile_TSS.log")
    shell:
        """
        plotProfile \
            --matrixFile {input.matrix} \
            --outFileName {output.profile} \
            --outFileNameData {output.data} \
            --perGroup \
            --refPointLabel TSS \
            --plotType se \
            --legendLocation upper-right \
            > {log} 2>&1
        """

rule computeMatrix_body:
    input:
        genes=os.path.join(REFERENCE_REGIONS_DIR, "genes.bed"),
        bws=expand(os.path.join(PER_SAMPLE_BIGWIG_DIR, "{sample}.dedup.bw"), sample=SAMPLES)
    output:
        matrix=os.path.join(GROUPED_QC_DIR, "matrix_all_bw_scalar.gz")
    params:
        body_length=config["deeptools"]["gene_body_length"],
        upstream=config["deeptools"]["gene_body_upstream"],
        downstream=config["deeptools"]["gene_body_downstream"]
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
        matrix=os.path.join(GROUPED_QC_DIR, "matrix_all_bw_scalar.gz")
    output:
        profile=os.path.join(GROUPED_QC_DIR, "all_samples_scalar_profile.pdf"),
        data=os.path.join(GROUPED_QC_DIR, "all_samples_scalar_profile.tab")
    log:
        os.path.join(LOG_DIR, "plotProfile_body.log")
    shell:
        """
        plotProfile \
            --matrixFile {input.matrix} \
            --outFileName {output.profile} \
            --outFileNameData {output.data} \
            --perGroup \
            --plotType se \
            --legendLocation upper-right \
            > {log} 2>&1
        """

# =====================================================
# 12b. Grouped Profile Plots (per qc_grouping)
# =====================================================
rule grouped_computeMatrix_TSS:
    input:
        tss=os.path.join(REFERENCE_REGIONS_DIR, "tss.bed"),
        bws=get_group_bigwigs
    output:
        matrix=os.path.join(GROUPED_QC_DIR, "{group}", "matrix_TSS.gz")
    params:
        labels=get_group_labels,
        upstream=config["deeptools"]["tss_upstream"],
        downstream=config["deeptools"]["tss_downstream"]
    threads: THREADS
    resources:
        mem_mb=config["resources"]["deeptools_mem_mb"]
    log:
        os.path.join(LOG_DIR, "{group}_computeMatrix_TSS.log")
    shell:
        """
        computeMatrix reference-point \
            --referencePoint TSS \
            -b {params.upstream} \
            -a {params.downstream} \
            -R {input.tss} \
            -S {input.bws} \
            --samplesLabel {params.labels} \
            --skipZeros \
            -o {output.matrix} \
            --binSize {BIN_SIZE} \
            --numberOfProcessors {threads} \
            > {log} 2>&1
        """

rule grouped_plotProfile_TSS:
    input:
        matrix=os.path.join(GROUPED_QC_DIR, "{group}", "matrix_TSS.gz")
    output:
        profile=os.path.join(GROUPED_QC_DIR, "{group}", "TSS_profile.pdf"),
        data=os.path.join(GROUPED_QC_DIR, "{group}", "TSS_profile.tab")
    log:
        os.path.join(LOG_DIR, "{group}_plotProfile_TSS.log")
    shell:
        """
        plotProfile \
            --matrixFile {input.matrix} \
            --outFileName {output.profile} \
            --outFileNameData {output.data} \
            --perGroup \
            --refPointLabel TSS \
            --plotType se \
            --legendLocation upper-right \
            -T "TSS Profile ({wildcards.group})" \
            > {log} 2>&1
        """

rule grouped_computeMatrix_body:
    input:
        genes=os.path.join(REFERENCE_REGIONS_DIR, "genes.bed"),
        bws=get_group_bigwigs
    output:
        matrix=os.path.join(GROUPED_QC_DIR, "{group}", "matrix_gene_body.gz")
    params:
        labels=get_group_labels,
        body_length=config["deeptools"]["gene_body_length"],
        upstream=config["deeptools"]["gene_body_upstream"],
        downstream=config["deeptools"]["gene_body_downstream"]
    threads: THREADS
    resources:
        mem_mb=config["resources"]["deeptools_mem_mb"]
    log:
        os.path.join(LOG_DIR, "{group}_computeMatrix_body.log")
    shell:
        """
        computeMatrix scale-regions \
            -R {input.genes} \
            -S {input.bws} \
            --samplesLabel {params.labels} \
            --regionBodyLength {params.body_length} \
            -b {params.upstream} \
            -a {params.downstream} \
            --binSize {BIN_SIZE} \
            --skipZeros \
            -o {output.matrix} \
            --numberOfProcessors {threads} \
            > {log} 2>&1
        """

rule grouped_plotProfile_body:
    input:
        matrix=os.path.join(GROUPED_QC_DIR, "{group}", "matrix_gene_body.gz")
    output:
        profile=os.path.join(GROUPED_QC_DIR, "{group}", "gene_body_profile.pdf"),
        data=os.path.join(GROUPED_QC_DIR, "{group}", "gene_body_profile.tab")
    log:
        os.path.join(LOG_DIR, "{group}_plotProfile_body.log")
    shell:
        """
        plotProfile \
            --matrixFile {input.matrix} \
            --outFileName {output.profile} \
            --outFileNameData {output.data} \
            --perGroup \
            --plotType se \
            --legendLocation upper-right \
            -T "Gene Body Profile ({wildcards.group})" \
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
            xls=os.path.join(PEAK_CALL_DIR, "{condition}", "{assay}", "{peak_type}", "{sample}_peaks.xls"),
            narrowPeak=os.path.join(PEAK_CALL_DIR, "{condition}", "{assay}", "{peak_type}", "{sample}_peaks.narrowPeak"),
            summits=os.path.join(PEAK_CALL_DIR, "{condition}", "{assay}", "{peak_type}", "{sample}_summits.bed"),
            broadPeak=os.path.join(PEAK_CALL_DIR, "{condition}", "{assay}", "{peak_type}", "{sample}_peaks.broadPeak"),
            gappedPeak=os.path.join(PEAK_CALL_DIR, "{condition}", "{assay}", "{peak_type}", "{sample}_peaks.gappedPeak")
        params:
            name="{sample}",
            outdir=lambda wc: os.path.join(PEAK_CALL_DIR, wc.condition, wc.assay, wc.peak_type),
            genome_size=config.get("macs3", {}).get("params", {}).get("genome_size", "hs"),
            pvalue=config.get("macs3", {}).get("params", {}).get("pvalue_cutoff", 0.01),
            keep_dup=config.get("macs3", {}).get("params", {}).get("keep_dup", "all"),
            broad_cutoff=config.get("macs3", {}).get("params", {}).get("broad_cutoff", 0.1),
            control_flag=lambda wc, input: f"-c {input.control}" if input.control else "",
            # Optional parameters
            format=config.get("macs3", {}).get("params", {}).get("format", "BAM"),
            bdg_flag="-B" if config.get("macs3", {}).get("params", {}).get("generate_bedgraph", False) else "",
            spmr_flag="--SPMR" if config.get("macs3", {}).get("params", {}).get("bedgraph_spmr", False) else "",
            model_flag="--nomodel --extsize {}".format(config.get("macs3", {}).get("params", {}).get("extsize", 200)) 
                       if config.get("macs3", {}).get("params", {}).get("no_model", False) else ""
        log:
            os.path.join(LOG_DIR, "{condition}_{assay}_{peak_type}_{sample}_macs3.log")
        wildcard_constraints:
            sample="|".join(MACS3_SAMPLES),
            peak_type="narrow|broad"
        shell:
            """
            mkdir -p {params.outdir}
            
            # Peak calling based on peak_type wildcard
            if [ "{wildcards.peak_type}" = "broad" ]; then
                # Broad peak calling
                macs3 callpeak \\
                    -t {input.treatment} \\
                    {params.control_flag} \\
                    -f {params.format} \\
                    -g {params.genome_size} \\
                    -n {params.name} \\
                    -p {params.pvalue} \\
                    --keep-dup {params.keep_dup} \\
                    --broad --broad-cutoff {params.broad_cutoff} \\
                    {params.bdg_flag} \\
                    {params.spmr_flag} \\
                    {params.model_flag} \\
                    --outdir {params.outdir} \\
                    > {log} 2>&1
                # Touch narrow files (not created in broad mode)
                touch {output.narrowPeak}
                touch {output.summits}
                
            else
                # Narrow peak calling
                macs3 callpeak \\
                    -t {input.treatment} \\
                    {params.control_flag} \\
                    -f {params.format} \\
                    -g {params.genome_size} \\
                    -n {params.name} \\
                    -p {params.pvalue} \\
                    --keep-dup {params.keep_dup} \\
                    {params.bdg_flag} \\
                    {params.spmr_flag} \\
                    {params.model_flag} \\
                    --outdir {params.outdir} \\
                    > {log} 2>&1
                # Touch broad files (not created in narrow mode)
                touch {output.broadPeak}
                touch {output.gappedPeak}
            fi
            """
# =====================================================
# 14. FRiP (Fraction of Reads in Peaks) Calculation
# =====================================================
if MACS3_SAMPLES:
    rule calculate_frip:
        input:
            bam=os.path.join(OUT_DEDUP, "{sample}.dedup.bam"),
            peaks=lambda wc: os.path.join(PEAK_CALL_DIR, wc.condition, wc.assay, wc.peak_type, 
                                         f"{wc.sample}_peaks.{'narrowPeak' if wc.peak_type == 'narrow' else 'broadPeak'}")
        output:
            frip=os.path.join(FRIP_DIR, "{condition}", "{assay}", "{peak_type}", "{sample}.frip.txt"),
            merged_peaks=os.path.join(FRIP_DIR, "{condition}", "{assay}", "{peak_type}", "{sample}.peaks.merged.bed")
        log:
            os.path.join(LOG_DIR, "{condition}_{assay}_{peak_type}_{sample}_frip.log")
        wildcard_constraints:
            sample="|".join(MACS3_SAMPLES),
            peak_type="narrow|broad"
        shell:
            """
            mkdir -p {FRIP_DIR}/{wildcards.condition}/{wildcards.assay}/{wildcards.peak_type}
            
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
            echo "Peak Type: {wildcards.peak_type}" >> {output.frip}
            echo "Total mapped reads : $TOTAL_READS" >> {output.frip}
            echo "Reads in peaks     : $READS_IN_PEAKS" >> {output.frip}
            echo "FRiP               : $FRIP" >> {output.frip}
            
            echo "FRiP calculation completed for {wildcards.sample} ({wildcards.peak_type})" > {log} 2>&1
            """
    rule aggregate_frip:
        input:
            frip_files=expand(os.path.join(FRIP_DIR, "{sample}.frip.txt"), sample=MACS3_SAMPLES)
        output:
            summary=os.path.join(FRIP_DIR, "all_samples_frip.tsv")
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
            
            # Log completion
            with open(log[0], 'w') as log_file:
                log_file.write("FRiP summary table created\n")
# =====================================================
# 15. IDR Analysis (Irreproducible Discovery Rate)
# =====================================================
if MACS3_SAMPLES_CONFIG and config.get("metadata_tsv"):
    rule idr_pairwise:
        input:
            peak1=lambda wc: os.path.join(PEAK_CALL_DIR, wc.condition, wc.assay, wc.peak_type,
                f"{get_sample_for_replicate_assay(wc.condition, wc.assay, wc.rep1)}_peaks."
                f"{'narrowPeak' if wc.peak_type == 'narrow' else 'broadPeak'}"),
            peak2=lambda wc: os.path.join(PEAK_CALL_DIR, wc.condition, wc.assay, wc.peak_type,
                f"{get_sample_for_replicate_assay(wc.condition, wc.assay, wc.rep2)}_peaks."
                f"{'narrowPeak' if wc.peak_type == 'narrow' else 'broadPeak'}")

        output:
            idr_output=os.path.join(IDR_DIR, "{condition}", "{assay}", "{peak_type}", "rep{rep1}_vs_rep{rep2}_idr.txt"),
            idr_plot=os.path.join(IDR_DIR, "{condition}", "{assay}", "{peak_type}", "rep{rep1}_vs_rep{rep2}_idr.txt.png")
        params:
            peak_file_type=lambda wc: "narrowPeak" if wc.peak_type == "narrow" else "broadPeak",
            rank_column=config.get("idr", {}).get("rank_column", "signal.value")
        log:
            os.path.join(LOG_DIR, "{condition}_{assay}_{peak_type}_rep{rep1}_vs{rep2}_idr.log")
        shell:
            """
            mkdir -p {IDR_DIR}/{wildcards.condition}/{wildcards.assay}/{wildcards.peak_type}
            
            idr --samples {input.peak1} {input.peak2} \\
                --input-file-type {params.peak_file_type} \\
                --rank {params.rank_column} \\
                --output-file {output.idr_output} \\
                --plot \\
                --log-output-file {log}
            """
    rule consensus_peaks:
        input:
            idr_files=lambda wc: [
                os.path.join(IDR_DIR, wc.condition, wc.assay, wc.peak_type, f"rep{r1}_vs_rep{r2}_idr.txt")
                for r1, r2 in get_replicate_pairs_for_assay(wc.condition, wc.assay)
            ],
            bam=lambda wc: os.path.join(OUT_DEDUP, f"{get_samples_for_condition_assay(wc.condition, wc.assay)[0]}.dedup.bam")
        output:
            bed=os.path.join(CONSENSUS_DIR, "{condition}", "{assay}", "{peak_type}", "consensus.bed"),
            frip=os.path.join(CONSENSUS_DIR, "{condition}", "{assay}", "{peak_type}", "consensus_frip.txt")
        log:
            os.path.join(LOG_DIR, "{condition}_{assay}_{peak_type}_consensus.log")
        shell:
            """
            mkdir -p {CONSENSUS_DIR}/{wildcards.condition}/{wildcards.assay}/{wildcards.peak_type}
            
            # Filter peaks with IDR >= 1.3 (more stringent threshold)
            # Merge all passing peaks from all pairwise comparisons
            cat {input.idr_files} | awk '$12 >= 1.3 {{print $1"\\t"$2"\\t"$3}}' | \\
              sort -k1,1 -k2,2n | \\
              bedtools merge -i stdin > {output.bed}
            
            # Calculate FRiP on consensus peaks
            TOTAL_READS=$(samtools view -c {input.bam})
            READS_IN_PEAKS=$(samtools view -b {input.bam} | \\
              bedtools intersect -u -a stdin -b {output.bed} | \\
              samtools view -c)
            FRIP=$(awk -v n="$READS_IN_PEAKS" -v d="$TOTAL_READS" 'BEGIN {{printf "%.5f", n/d}}')
            
            echo "Condition: {wildcards.condition}" > {output.frip}
            echo "Assay: {wildcards.assay}" >> {output.frip}
            echo "Peak Type: {wildcards.peak_type}" >> {output.frip}
            echo "Total reads: $TOTAL_READS" >> {output.frip}
            echo "Reads in consensus peaks: $READS_IN_PEAKS" >> {output.frip}
            echo "FRiP (consensus): $FRIP" >> {output.frip}
            echo "Number of consensus peaks: $(wc -l < {output.bed})" >> {output.frip}
            
            echo "Consensus peaks generated for {wildcards.condition} {wildcards.assay} {wildcards.peak_type}" > {log} 2>&1
            """



# =====================================================
# 16. HOMER Motif Analysis (Optional)
# =====================================================
if RUN_HOMER and MACS3_SAMPLES_CONFIG and config.get("metadata_tsv"):
    rule homer_motif_analysis:
        input:
            consensus=os.path.join(CONSENSUS_DIR, "{condition}", "{assay}", "{peak_type}", "consensus.bed"),
            genome=GENOME_FASTA
        output:
            html=os.path.join(HOMER_DIR, "{condition}", "{assay}", "{peak_type}", "homerResults.html"),
            motifs=directory(os.path.join(HOMER_DIR, "{condition}", "{assay}", "{peak_type}"))

        params:
            size=config.get("homer_params", {}).get("size", 200),
            length=config.get("homer_params", {}).get("length", "8,10,12"),
            mask="-mask" if config.get("homer_params", {}).get("mask", True) else ""
        threads: THREADS
        log:
            os.path.join(LOG_DIR, "{condition}_{assay}_{peak_type}_homer.log")
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
                os.path.join(PER_SAMPLE_BIGWIG_DIR, "{sample}.dedup.bw"),
                sample=get_samples_for_condition_assay(wc.condition, wc.assay)
            )
        output:
            mean_bw=os.path.join(AVERAGED_BIGWIG_DIR, "{condition}", "{assay}_mean.bw")

        threads: THREADS
        log:
            os.path.join(LOG_DIR, "{condition}_{assay}_bigwigAverage.log")
        shell:
            """
            mkdir -p {AVERAGED_BIGWIG_DIR}/{wildcards.condition}
            
            bigwigAverage \
                -b {input.bigwigs} \
                -o {output.mean_bw} \
                -p {threads} \
                > {log} 2>&1
            """

# =====================================================
# 18. Normalize BigWigs (ChIP/Input Log2 Ratio)
# =====================================================
if MACS3_SAMPLES_CONFIG and config.get("metadata_tsv"):
    rule calculate_normalized_bigwig:
        input:
            chip_bw=os.path.join(AVERAGED_BIGWIG_DIR, "{condition}", "{assay}_mean.bw"),
            input_bw=os.path.join(AVERAGED_BIGWIG_DIR, "{condition}", "input_mean.bw")
        output:
            normalized_bw=os.path.join(NORMALIZED_BIGWIG_DIR, "{condition}", "{assay}_log2_over_input.bw")

        params:
            pseudocount=1
        threads: THREADS
        log:
            os.path.join(LOG_DIR, "{condition}_{assay}_bigwigCompare.log")
        shell:
            """
            mkdir -p {NORMALIZED_BIGWIG_DIR}/{wildcards.condition}

            
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
            tss=os.path.join(REFERENCE_REGIONS_DIR, "tss.bed"),
            bw=os.path.join(NORMALIZED_BIGWIG_DIR, "{condition}", "{assay}_log2_over_input.bw")
        output:
            matrix=os.path.join(ENRICHMENT_PROFILES_DIR, "{condition}", "{assay}", "TSS_log2.mat.gz"),
            profile=report(
                os.path.join(ENRICHMENT_PROFILES_DIR, "{condition}", "{assay}", "TSS_log2_profile.pdf"),
                caption="report/tss_profile.rst",
                category="Profile Plots",
                labels={"condition": "{condition}", "assay": "{assay}", "type": "TSS_normalized"}
            ),
            data=os.path.join(ENRICHMENT_PROFILES_DIR, "{condition}", "{assay}", "TSS_log2_profile.tab")

        params:
            upstream=config["deeptools"]["tss_upstream"],
            downstream=config["deeptools"]["tss_downstream"]
        threads: THREADS
        log:
            os.path.join(LOG_DIR, "{condition}_{assay}_TSS_log2_profile.log")
        shell:
            """
            mkdir -p {ENRICHMENT_PROFILES_DIR}/{wildcards.condition}/{wildcards.assay}

            
            computeMatrix reference-point \
                --referencePoint TSS \
                -b {params.upstream} -a {params.downstream} \
                --binSize {BIN_SIZE} \
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
                -out {output.profile} \
                --outFileNameData {output.data} \
                >> {log} 2>&1
            """
    # Gene body profiles - normalized (log2)
    rule profile_genes_normalized:
        input:
            genes=os.path.join(REFERENCE_REGIONS_DIR, "genes.bed"),
            bw=os.path.join(NORMALIZED_BIGWIG_DIR, "{condition}", "{assay}_log2_over_input.bw")
        output:
            matrix=os.path.join(ENRICHMENT_PROFILES_DIR, "{condition}", "{assay}", "genes_log2.mat.gz"),
            profile=os.path.join(ENRICHMENT_PROFILES_DIR, "{condition}", "{assay}", "genes_log2_profile.pdf"),
            data=os.path.join(ENRICHMENT_PROFILES_DIR, "{condition}", "{assay}", "genes_log2_profile.tab")

        params:
            upstream=config["deeptools"]["gene_body_upstream"],
            downstream=config["deeptools"]["gene_body_downstream"],
            body_length=config["deeptools"]["gene_body_length"]
        threads: THREADS
        log:
            os.path.join(LOG_DIR, "{condition}_{assay}_genes_log2_profile.log")
        shell:
            """
            mkdir -p {ENRICHMENT_PROFILES_DIR}/{wildcards.condition}/{wildcards.assay}

            
            computeMatrix scale-regions \
                -b {params.upstream} -a {params.downstream} \
                --regionBodyLength {params.body_length} \
                --binSize {BIN_SIZE} \
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
                -out {output.profile} \
                --outFileNameData {output.data} \
                >> {log} 2>&1
            """
    # TSS profiles - averaged (mean signal)
    rule profile_tss_mean:
        input:
            tss=os.path.join(REFERENCE_REGIONS_DIR, "tss.bed"),
            bw=os.path.join(AVERAGED_BIGWIG_DIR, "{condition}", "{assay}_mean.bw")
        output:
            matrix=os.path.join(ENRICHMENT_PROFILES_DIR, "{condition}", "{assay}", "TSS_mean.mat.gz"),
            profile=os.path.join(ENRICHMENT_PROFILES_DIR, "{condition}", "{assay}", "TSS_mean_profile.pdf"),
            data=os.path.join(ENRICHMENT_PROFILES_DIR, "{condition}", "{assay}", "TSS_mean_profile.tab")

        params:
            upstream=config["deeptools"]["tss_upstream"],
            downstream=config["deeptools"]["tss_downstream"]
        threads: THREADS
        log:
            os.path.join(LOG_DIR, "{condition}_{assay}_TSS_mean_profile.log")
        shell:
            """
            mkdir -p {ENRICHMENT_PROFILES_DIR}/{wildcards.condition}/{wildcards.assay}

            
            computeMatrix reference-point \
                --referencePoint TSS \
                -b {params.upstream} -a {params.downstream} \
                --binSize {BIN_SIZE} \
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
                -out {output.profile} \
                --outFileNameData {output.data} \
                >> {log} 2>&1
            """
    # Gene body profiles - averaged (mean signal)
    rule profile_genes_mean:
        input:
            genes=os.path.join(REFERENCE_REGIONS_DIR, "genes.bed"),
            bw=os.path.join(AVERAGED_BIGWIG_DIR, "{condition}", "{assay}_mean.bw")
        output:
            matrix=os.path.join(ENRICHMENT_PROFILES_DIR, "{condition}", "{assay}", "genes_mean.mat.gz"),
            profile=os.path.join(ENRICHMENT_PROFILES_DIR, "{condition}", "{assay}", "genes_mean_profile.pdf"),
            data=os.path.join(ENRICHMENT_PROFILES_DIR, "{condition}", "{assay}", "genes_mean_profile.tab")

        params:
            upstream=config["deeptools"]["gene_body_upstream"],
            downstream=config["deeptools"]["gene_body_downstream"],
            body_length=config["deeptools"]["gene_body_length"]
        threads: THREADS
        log:
            os.path.join(LOG_DIR, "{condition}_{assay}_genes_mean_profile.log")
        shell:
            """
            mkdir -p {ENRICHMENT_PROFILES_DIR}/{wildcards.condition}/{wildcards.assay}

            
            computeMatrix scale-regions \
                -b {params.upstream} -a {params.downstream} \
                --regionBodyLength {params.body_length} \
                --binSize {BIN_SIZE} \
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
                -out {output.profile} \
                --outFileNameData {output.data} \
                >> {log} 2>&1
            """
# =====================================================
# 20. Generate Comprehensive Excel Summary
# =====================================================
rule generate_analysis_summary:
    input:
        bigwigs=expand(os.path.join(PER_SAMPLE_BIGWIG_DIR, "{sample}.dedup.bw"), sample=SAMPLES),
        frip=ALL_FRIP if ALL_FRIP else [],
        idr=ALL_IDR if ALL_IDR else [],
        consensus=ALL_CONSENSUS_PEAKS if ALL_CONSENSUS_PEAKS else [],
        qc=ALL_DEEPTOOLS_QC,
        metadata=config.get("metadata_tsv", "samples_metadata_example.tsv")
    output:
        summary=ANALYSIS_SUMMARY_PATH
    log:
        os.path.join(LOG_DIR, "generate_summary.log")
    run:
        import pandas as pd
        from datetime import datetime
        import glob
        
        print("Generating comprehensive analysis summary...", file=open(log[0], 'w'))
        
        # Load metadata
        metadata_df = pd.read_csv(input.metadata, sep="\t") if os.path.exists(input.metadata) else None
        
        # Pre-calculate groups for consistency
        condition_assay_all = {}  # All assays (ChIP + Input)
        condition_assay_chip = {} # Only ChIP assays
        for sample, info in MACS3_SAMPLES_CONFIG.items():
            cond = info.get("condition", "default")
            assay = info.get("assay", "default")
            key = (cond, assay)
            if key not in condition_assay_all:
                condition_assay_all[key] = []
            condition_assay_all[key].append(sample)
            
            if info.get("mark_type"):
                if key not in condition_assay_chip:
                    condition_assay_chip[key] = []
                condition_assay_chip[key].append(sample)
        
        with pd.ExcelWriter(output.summary, engine='openpyxl') as writer:
            # Sheet 1: Raw FASTQ Files
            raw_fastq_data = []
            for sample in SAMPLES:
                if SEQ_TYPE == "paired":
                    raw_fastq_data.append({
                        "Sample": sample,
                        "R1_FASTQ": f"{RAW_READS_DIR}/{sample}_R1.fastq.gz",
                        "R2_FASTQ": f"{RAW_READS_DIR}/{sample}_R2.fastq.gz"
                    })
                else:
                    raw_fastq_data.append({
                        "Sample": sample,
                        "FASTQ": f"{RAW_READS_DIR}/{sample}.fastq.gz"
                    })
            pd.DataFrame(raw_fastq_data).to_excel(writer, sheet_name="Raw_FASTQ", index=False)
            
            # Sheet 2: Trimmed FASTQ Files
            trimmed_fastq_data = []
            for sample in SAMPLES:
                if SEQ_TYPE == "paired":
                    trimmed_fastq_data.append({
                        "Sample": sample,
                        "R1_Trimmed": f"{TRIM_DIR}/{sample}_R1_val_1.fq.gz",
                        "R2_Trimmed": f"{TRIM_DIR}/{sample}_R2_val_2.fq.gz",
                        "R1_Report": f"{TRIM_DIR}/{sample}_R1.fastq.gz_trimming_report.txt",
                        "R2_Report": f"{TRIM_DIR}/{sample}_R2.fastq.gz_trimming_report.txt"
                    })
                else:
                    trimmed_fastq_data.append({
                        "Sample": sample,
                        "Trimmed": f"{TRIM_DIR}/{sample}_trimmed.fq.gz",
                        "Report": f"{TRIM_DIR}/{sample}.fastq.gz_trimming_report.txt"
                    })
            pd.DataFrame(trimmed_fastq_data).to_excel(writer, sheet_name="Trimmed_FASTQ", index=False)
            
            # Sheet 3: BAM Files
            bam_data = []
            for sample in SAMPLES:
                bam_data.append({
                    "Sample": sample,
                    "Sorted_BAM": f"{OUT_ALIGN}/{sample}.sorted.bam",
                    "Sorted_BAI": f"{OUT_ALIGN}/{sample}.sorted.bam.bai",
                    "Dedup_BAM": f"{OUT_DEDUP}/{sample}.dedup.bam",
                    "Dedup_BAI": f"{OUT_DEDUP}/{sample}.dedup.bam.bai",
                    "BigWig": os.path.join(PER_SAMPLE_BIGWIG_DIR, f"{sample}.dedup.bw")
                })
            pd.DataFrame(bam_data).to_excel(writer, sheet_name="BAM_Files", index=False)
            
            # Sheet 4: QC Metrics
            qc_metrics_data = []
            for sample in SAMPLES:
                qc_metrics_data.append({
                    "Sample": sample,
                    "Flagstat": f"{OUT_QC}/{sample}.flagstat.txt",
                    "Stats": f"{OUT_QC}/{sample}.stats.txt",
                    "PBC": f"{CHIP_QC_DIR}/{sample}.pbc.txt"
                })
            pd.DataFrame(qc_metrics_data).to_excel(writer, sheet_name="QC_Metrics", index=False)
            
            # Sheet 5: SPP Analysis (if available)
            if RUN_SPP:
                spp_data = []
                for sample in SAMPLES:
                    spp_data.append({
                        "Sample": sample,
                        "SPP_QC": f"{CHIP_QC_DIR}/{sample}_spp.qc.txt",
                        "AVP_Plot": f"{CHIP_QC_DIR}/{sample}_avp.pdf"
                    })
                pd.DataFrame(spp_data).to_excel(writer, sheet_name="SPP_Analysis", index=False)
            
            # Sheet 6: MACS3 Peak Calling
            if MACS3_SAMPLES_CONFIG:
                macs3_data = []
                for sample, info in MACS3_SAMPLES_CONFIG.items():
                    if info.get("mark_type"):
                        condition = info.get("condition", "default")
                        assay = info.get("assay", "default")
                        peak_types = get_peak_types_for_sample(sample)
                        for peak_type in peak_types:
                            peak_ext = "narrowPeak" if peak_type == "narrow" else "broadPeak"
                            macs3_data.append({
                                "Sample": sample,
                                "Condition": condition,
                                "Replicate": info.get("replicate", ""),
                                "Assay": assay,
                                "Peak_Type": peak_type,
                                "Control": info.get("control", ""),
                                "Treatment_BAM": f"{OUT_DEDUP}/{sample}.dedup.bam",
                                "Peaks_File": f"macs3_results/{condition}/{assay}/{peak_type}/{sample}_peaks.{peak_ext}",
                                "Summits_or_Gapped": f"macs3_results/{condition}/{assay}/{peak_type}/{sample}_{'summits.bed' if peak_type=='narrow' else 'peaks.gappedPeak'}",
                                "XLS_Report": f"macs3_results/{condition}/{assay}/{peak_type}/{sample}_peaks.xls"
                            })
                if macs3_data:
                    pd.DataFrame(macs3_data).to_excel(writer, sheet_name="MACS3_Peaks", index=False)

            
            # Sheet 7: FRiP Analysis
            if MACS3_SAMPLES:
                frip_data = []
                for sample in MACS3_SAMPLES:
                    info = MACS3_SAMPLES_CONFIG[sample]
                    condition = info.get("condition", "default")
                    assay = info.get("assay", "default")
                    peak_types = get_peak_types_for_sample(sample)
                    
                    for peak_type in peak_types:
                        frip_file = os.path.join(FRIP_DIR, condition, assay, peak_type, f"{sample}.frip.txt")
                        if os.path.exists(frip_file):
                            with open(frip_file) as f:
                                content = f.read()
                                try:
                                    total = content.split("Total mapped reads")[1].split(":")[1].strip().split()[0]
                                    in_peaks = content.split("Reads in peaks")[1].split(":")[1].strip().split()[0]
                                    frip = content.split("FRiP")[1].split(":")[1].strip().split()[0]
                                except:
                                    total = in_peaks = frip = "N/A"
                            
                            frip_data.append({
                                "Sample": sample,
                                "Assay": assay,
                                "Peak_Type": peak_type,
                                "Total_Reads": total,
                                "Reads_in_Peaks": in_peaks,
                                "FRiP": frip,
                                "FRiP_File": frip_file
                            })
                if frip_data:
                    pd.DataFrame(frip_data).to_excel(writer, sheet_name="FRiP_Analysis", index=False)

            
            # Sheet 8: IDR Comparisons
            if MACS3_SAMPLES_CONFIG:
                idr_data = []
                for (condition, assay), samples in condition_assay_chip.items():
                    peak_types = get_all_peak_types_for_condition_assay(condition, assay)
                    if len(samples) >= 2:
                        for sample1, sample2 in combinations(sorted(samples), 2):
                            rep1 = MACS3_SAMPLES_CONFIG[sample1]["replicate"]
                            rep2 = MACS3_SAMPLES_CONFIG[sample2]["replicate"]
                            for peak_type in peak_types:
                                idr_file = os.path.join(IDR_DIR, condition, assay, peak_type, f"rep{rep1}_vs_rep{rep2}_idr.txt")
                                idr_data.append({
                                    "Condition": condition,
                                    "Assay": assay,
                                    "Peak_Type": peak_type,
                                    "Rep1": rep1,
                                    "Rep2": rep2,
                                    "IDR_File": idr_file,
                                    "IDR_Plot": f"{idr_file}.png"
                                })
                if idr_data:
                    pd.DataFrame(idr_data).to_excel(writer, sheet_name="IDR_Comparisons", index=False)

            
            # Sheet 9: Consensus Peaks
            if ALL_CONSENSUS_PEAKS:
                consensus_data = []
                for (condition, assay), samples in condition_assay_chip.items():
                    peak_types = get_all_peak_types_for_condition_assay(condition, assay)
                    for peak_type in peak_types:
                        consensus_data.append({
                            "Condition": condition,
                            "Assay": assay,
                            "Peak_Type": peak_type,
                            "Consensus_BED": os.path.join(CONSENSUS_DIR, condition, assay, peak_type, "consensus.bed"),
                            "FRiP_Report": os.path.join(CONSENSUS_DIR, condition, assay, peak_type, "consensus_frip.txt")
                        })
                if consensus_data:
                    pd.DataFrame(consensus_data).to_excel(writer, sheet_name="Consensus_Peaks", index=False)
            
            # Sheet 10: HOMER Motifs (if available)
            if RUN_HOMER and ALL_HOMER:
                homer_data = []
                for (condition, assay), samples in condition_assay_chip.items():
                    peak_types = get_all_peak_types_for_condition_assay(condition, assay)
                    for peak_type in peak_types:
                        homer_data.append({
                            "Condition": condition,
                            "Assay": assay,
                            "Peak_Type": peak_type,
                            "HOMER_Dir": os.path.join(HOMER_DIR, condition, assay, peak_type),
                            "HTML_Report": os.path.join(HOMER_DIR, condition, assay, peak_type, "homerResults.html")
                        })
                if homer_data:
                    pd.DataFrame(homer_data).to_excel(writer, sheet_name="HOMER_Motifs", index=False)
            
            # Sheet 11: BigWig Files
            if MACS3_SAMPLES_CONFIG:
                bigwig_data = []
                for (condition, assay), samples in condition_assay_all.items():
                    is_chip = any(MACS3_SAMPLES_CONFIG[s].get("mark_type") for s in samples)
                    bigwig_data.append({
                        "Condition": condition,
                        "Assay": assay,
                        "Averaged_BW": os.path.join(AVERAGED_BIGWIG_DIR, condition, f"{assay}_mean.bw"),
                        "Normalized_BW": os.path.join(NORMALIZED_BIGWIG_DIR, condition, f"{assay}_log2_over_input.bw") if is_chip else "N/A"
                    })
                if bigwig_data:
                    pd.DataFrame(bigwig_data).to_excel(writer, sheet_name="BigWig_Files", index=False)

            
            # Sheet 12: Profile Plots
            if ALL_CONDITION_PROFILES:
                profile_data = []
                for (condition, assay), samples in condition_assay_chip.items():
                    profile_data.append({
                        "Condition": condition,
                        "Assay": assay,
                        "TSS_Normalized_Plot": os.path.join(ENRICHMENT_PROFILES_DIR, condition, assay, "TSS_log2_profile.pdf"),
                        "TSS_Mean_Plot": os.path.join(ENRICHMENT_PROFILES_DIR, condition, assay, "TSS_mean_profile.pdf"),
                        "Genes_Normalized_Plot": os.path.join(ENRICHMENT_PROFILES_DIR, condition, assay, "genes_log2_profile.pdf"),
                        "Genes_Mean_Plot": os.path.join(ENRICHMENT_PROFILES_DIR, condition, assay, "genes_mean_profile.pdf")
                    })
                if profile_data:
                    pd.DataFrame(profile_data).to_excel(writer, sheet_name="Profile_Plots", index=False)

            
            # Sheet 13: QC Outputs
            qc_data = []
            for grp in DEEPTOOLS_GROUPS:
                qc_data.extend([
                    {"Analysis": f"Correlation Heatmap ({grp})", "File": os.path.join(GROUPED_QC_DIR, grp, "correlation_heatmap.pdf")},
                    {"Analysis": f"Correlation Scatter ({grp})", "File": os.path.join(GROUPED_QC_DIR, grp, "correlation_scatterplot.pdf")},
                    {"Analysis": f"PCA Plot ({grp})", "File": os.path.join(GROUPED_QC_DIR, grp, "pca.pdf")},
                    {"Analysis": f"Fingerprint Plot ({grp})", "File": os.path.join(GROUPED_QC_DIR, grp, "fingerprints.pdf")},
                    {"Analysis": f"Coverage Histogram ({grp})", "File": os.path.join(GROUPED_QC_DIR, grp, "coverage_histogram.pdf")},
                ])
            qc_data.append({"Analysis": "MultiQC Report", "File": os.path.join(OUT_QC, "alignment_qc_report.html")})
            
            pd.DataFrame(qc_data).to_excel(writer, sheet_name="QC_Outputs", index=False)
        
        print(f"✅ Analysis summary saved to {output.summary}", file=open(log[0], 'a'))

# =====================================================
# 21. Generate Automated MkDocs Results Summary
# =====================================================
rule generate_mkdocs_results:
    input:
        frip=ALL_FRIP if ALL_FRIP else [],
        idr=ALL_IDR if ALL_IDR else [],
        consensus=ALL_CONSENSUS_PEAKS if ALL_CONSENSUS_PEAKS else [],
        qc=ALL_DEEPTOOLS_QC,
        metadata=config.get("metadata_tsv", "samples_metadata_example.tsv")
    output:
        md="docs/results/latest_run.md"
    log:
        os.path.join(LOG_DIR, "generate_mkdocs_results.log")
    run:
        import pandas as pd
        import os
        from datetime import datetime

        with open(output.md, 'w') as f:
            f.write(f"# Latest Analysis Results\n\n")
            f.write(f"*Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}*\n\n")
            
            # --- Project Summary ---
            f.write("## 📋 Project Summary\n\n")
            f.write(f"- **Sequencing Type**: {SEQ_TYPE.upper()}\n")
            f.write(f"- **Samples Processed**: {len(SAMPLES)}\n")
            f.write(f"- **Genome**: {config['genome'].get('fasta', 'N/A')}\n\n")

            # --- Sample Groups (Mermaid) ---
            f.write("## 🧬 Sample Groups\n\n")
            f.write("```mermaid\ngraph LR\n")
            # Build groupings
            cond_assay_all = {}
            cond_assay_chip = {}
            for s, info in MACS3_SAMPLES_CONFIG.items():
                c = info.get('condition', 'default')
                a = info.get('assay', 'default')
                k = (c, a)
                if k not in cond_assay_all: cond_assay_all[k] = []
                cond_assay_all[k].append(s)
                if info.get('mark_type'):
                    if k not in cond_assay_chip: cond_assay_chip[k] = []
                    cond_assay_chip[k].append(s)

            for (cond, assay), samples in cond_assay_chip.items():
                f.write(f'    subgraph \"{cond}_{assay}\"\n')
                for s in samples:
                    f.write(f'        {s}\n')
                f.write('    end\n')
            f.write("```\n\n")

            # --- FRiP Analysis ---
            if ALL_FRIP:
                f.write("## 📏 Peak Enrichment (FRiP)\n\n")
                f.write("| Sample | Assay | Type | FRiP | Total Reads |\n")
                f.write("| :--- | :--- | :--- | :--- | :--- |\n")
                for sample in MACS3_SAMPLES:
                    info = MACS3_SAMPLES_CONFIG[sample]
                    cond = info.get("condition", "default")
                    assay = info.get("assay", "default")
                    # Use the helper function if available, else assume narrow/broad
                    try:
                        ptypes = get_peak_types_for_sample(sample)
                    except:
                        ptypes = ["narrow"]
                    for ptype in ptypes:
                        frip_file = os.path.join(FRIP_DIR, cond, assay, ptype, f"{sample}.frip.txt")
                        if os.path.exists(frip_file):
                            with open(frip_file) as rf:
                                content = rf.read()
                                try:
                                    total = content.split("Total mapped reads")[1].split(":")[1].strip().split()[0]
                                    frip = content.split("FRiP")[1].split(":")[1].strip().split()[0]
                                    f.write(f"| {sample} | {assay} | {ptype} | **{frip}** | {total} |\n")
                                except:
                                    pass
                f.write("\n")

            # --- Visual Tracks ---
            f.write("## 🛰️ Visualization Tracks\n\n")
            f.write("!!! info \"Technical Note\"\n")
            f.write("    The averaged tracks provide a non-biased view of enrichment across replicates.\n\n")
            f.write("| Condition | Assay | Averaged Track (BW) | Normalized Track (log2) |\n")
            f.write("| :--- | :--- | :--- | :--- |\n")
            for (cond, assay), samples in cond_assay_all.items():
                is_chip = any(MACS3_SAMPLES_CONFIG[s].get("mark_type") for s in samples)
                avg_bw = f"../../03_deeptools_output/3.3_averaged_and_input_normalized_enrichment_outputs/3.3.1_replicate_averaged_bigwigs/{cond}/{assay}_mean.bw"
                norm_bw = f"../../03_deeptools_output/3.3_averaged_and_input_normalized_enrichment_outputs/3.3.2_chip_over_input_log2_bigwigs/{cond}/{assay}_log2_over_input.bw" if is_chip else "N/A"
                f.write(f"| {cond} | {assay} | [View Avg]({avg_bw}) | {'[View Norm]('+norm_bw+')' if is_chip else 'N/A'} |\n")
            f.write("\n")

            # --- Quality Control Plots ---
            f.write("## 📊 DeepTools QC Plots\n\n")
            f.write("QC plots are generated for config-defined groups (all samples, per-condition, etc.).\n\n")
            f.write("| Group | Metric | Static Output |\n")
            f.write("| :--- | :--- | :--- |\n")
            for grp in DEEPTOOLS_GROUPS:
                f.write(f"| **{grp}** | Correlation Heatmap | [View PDF](../../03_deeptools_output/3.2_replicatewise_qc_profiles/{grp}/correlation_heatmap.pdf) |\n")
                f.write(f"| | PCA Plot | [View PDF](../../03_deeptools_output/3.2_replicatewise_qc_profiles/{grp}/pca.pdf) |\n")
                f.write(f"| | Fingerprints | [View PDF](../../03_deeptools_output/3.2_replicatewise_qc_profiles/{grp}/fingerprints.pdf) |\n")
                f.write(f"| | Coverage Histogram | [View PDF](../../03_deeptools_output/3.2_replicatewise_qc_profiles/{grp}/coverage_histogram.pdf) |\n")
            f.write("\n")

            # --- Profile Plots ---
            f.write("## 📈 Enrichment Profiles\n\n")
            for (cond, assay), samples in cond_assay_chip.items():
                f.write(f"### {cond} - {assay}\n")
                f.write(f"**TSS Profile**: [Static PDF](../../03_deeptools_output/3.3_averaged_and_input_normalized_enrichment_outputs/3.3.3_tss_gene_body_enrichment_profiles/{cond}/{assay}/TSS_log2_profile.pdf)\n")
                f.write(f"- Mean Signal: [Static PDF](../../03_deeptools_output/3.3_averaged_and_input_normalized_enrichment_outputs/3.3.3_tss_gene_body_enrichment_profiles/{cond}/{assay}/TSS_mean_profile.pdf)\n\n")
                
                f.write(f"**Gene Body Profile**: [Static PDF](../../03_deeptools_output/3.3_averaged_and_input_normalized_enrichment_outputs/3.3.3_tss_gene_body_enrichment_profiles/{cond}/{assay}/genes_log2_profile.pdf)\n")
                f.write(f"- Mean Signal: [Static PDF](../../03_deeptools_output/3.3_averaged_and_input_normalized_enrichment_outputs/3.3.3_tss_gene_body_enrichment_profiles/{cond}/{assay}/genes_mean_profile.pdf)\n\n")
