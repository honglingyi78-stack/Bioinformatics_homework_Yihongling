# Week 5 Interpretation — DESeq2 Differential Expression

**Comparison:** treated versus control with batch in the design (`~ batch + condition`); 989 of 1,000 genes passed pre-filtering (≥10 counts in ≥3 samples).

**Strongest QC observation:** The VST-based PCA separates samples by condition along PC1 (24% of variance) while the three batches intermingle — a balanced design with no outlier sample.

**Significant genes:** At `padj < 0.05` and `|log2FoldChange| ≥ 1`, 60 genes were significant: 36 up-regulated and 24 down-regulated in treated samples.

**Biological interpretation:** More genes are induced than repressed, suggesting the treatment activates a coordinated transcriptional program; the largest effect (Gene0008 ≈ 3.8-fold) marks strong responders.

**Limitation:** Gene identifiers are anonymized, so no functional annotation is possible, and two replicates per batch–condition cell limit power for small effects.

**AI use:** AI drafted the analysis script and this summary; I independently verified sample identity, coefficient direction, and significance thresholds against the raw data and result table.
