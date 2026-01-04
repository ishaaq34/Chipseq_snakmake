Gene Body Coverage
==================

ChIP signal across gene bodies (TSS to TES) and flanking regions.

The plot shows:

* **Upstream**: 3kb before TSS
* **Gene body**: Scaled to 5kb
* **Downstream**: 3kb after TES
* **Y-axis**: Log2(ChIP/Input) normalized signal

Expected patterns:

* **Promoter marks**: Peak at TSS, decline in gene body
* **Elongation marks** (H3K36me3): Enrichment throughout gene body
* **Repressive marks**: Variable patterns depending on genomic context

Generated using DeepTools computeMatrix (scale-regions mode) and plotProfile.
