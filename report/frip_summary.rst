FRiP Summary Table
==================

Fraction of Reads in Peaks (FRiP) for all samples.

FRiP represents the proportion of reads falling within called peaks and is a key quality metric:

* FRiP > 0.01 (1%): Acceptable for broad marks
* FRiP > 0.05 (5%): Good for most marks
* FRiP > 0.10 (10%): Excellent signal enrichment

Higher values indicate better IP quality and enrichment.

Calculated using bedtools intersect on deduplicated BAM files.
