# Week 5 — AI Verification Log

## Task given to the AI

> "Based on the Week 5 starter script, produce a complete, runnable bulk
> RNA-seq differential-expression analysis: import the count matrix and
> metadata, validate the inputs, build a DESeq2 object with
> `design = ~ batch + condition`, pre-filter, run DESeq(), extract the
> treated-vs-control contrast, apply apeglm shrinkage, generate a PCA and a
> DE plot, export the full results table, save the fitted object and
> `sessionInfo()`, and write a 100–150 word interpretation."

## What AI generated

- `week5_deseq2_analysis.R` — the complete analysis script (adapted from the
  starter; TODOs completed).
- This log and the draft interpretation in `week5_interpretation.md`.

## What was verified

| Check | How it was verified | Result |
|---|---|---|
| Counts are non-negative integers | `stopifnot` in script: `all(counts >= 0)`, `counts == round(counts)` | Passed |
| Count-matrix columns = metadata row names, same order | `identical(colnames(counts), rownames(coldata))` | Passed (12/12) |
| No duplicated sample IDs | `!anyDuplicated(rownames(coldata))` | Passed |
| `condition` / `batch` levels | `table(batch, condition)` printed; levels checked | A/B/C × control/treated, balanced 2×2 |
| `control` is reference level | `levels(coldata$condition)[1] == "control"` | Passed |
| Design formula | `design(dds)` printed `~ batch + condition` | Confirmed |
| Coefficient name | `resultsNames(dds)` inspected → `condition_treated_vs_control`; script stops if absent | Confirmed |
| Contrast direction | Contrast `c("condition","treated","control")`; positive LFC = up in treated | Confirmed from top rows |
| Thresholds | `padj < 0.05` (alpha = 0.05 in `results()`) and `abs(log2FoldChange) >= 1` | Used for `significant` flag |
| Result table | 989 rows (all genes, incl. nonsignificant), 8 columns; 60 significant (36 up / 24 down) | Verified by re-reading CSV |
| Saved object | `readRDS()` returns `DESeqDataSet`, 989 × 12, design intact | Verified |
| Plots | `week5_pca.png`, `week5_de_plot.png` inspected after rendering | Axes, labels, thresholds OK |
| `sessionInfo()` | Saved to `session_info.txt` (R 4.5.1, DESeq2 1.48.1, apeglm 1.30.0) | Saved |

## AI-generated error and revision (documented)

- **Error:** `vst(dds, blind = FALSE)` failed with
  `less than 'nsub' rows ... use varianceStabilizingTransformation directly`.
  Cause: after pre-filtering, 989 genes remain, below the default `nsub = 1000`
  required by the fast `vst()` approximation.
- **Fix:** replaced with `varianceStabilizingTransformation(dds, blind = FALSE)`,
  which computes the same VST without the subset approximation. The script was
  re-run end-to-end and the PCA regenerated from the corrected transformation.

## Boundary of AI involvement

- The AI wrote and executed the code; the statistical model, contrast, and
  thresholds were taken from the assignment instructions.
- Sample identity (metadata ↔ count columns), coefficient direction
  (treated up = positive LFC), and the significance counts were independently
  re-checked against the raw inputs and exported table, not taken on trust.
- The biological interpretation is a draft; final claims and wording remain
  the student's responsibility.
