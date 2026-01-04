TSS Enrichment Profile
=======================

ChIP signal enrichment around transcription start sites (TSS).

The plot shows:

* **X-axis**: Distance from TSS (-3kb to +3kb)
* **Y-axis**: Log2(ChIP/Input) normalized signal

Expected patterns:

* **Active marks** (H3K4me3, H3K9ac, H3K27ac): Sharp peak at TSS
* **Repressive marks** (H3K27me3, H3K9me3): Depletion at TSS, enrichment in gene body
* **Enhancer marks** (H3K4me1): Enrichment at distal sites

Generated using DeepTools computeMatrix and plotProfile on averaged replicates.
