# ================================================================
# Week 5 Homework — Bulk RNA-seq differential expression with DESeq2
# From Verified Counts to an Interpretable DESeq2 Result
#
# This script is adapted from Week5_Homework_Starter.R and completes
# every step of the required workflow:
#   1. Import counts + metadata
#   2. Verify data integrity (non-negative integers, column/row match,
#      factor levels, control as reference)
#   3. DESeq2 object with design = ~ batch + condition
#   4. Pre-filter (>= 10 counts in >= 3 samples)
#   5. DESeq() + inspect resultsNames(dds)
#   6. Extract treated-vs-control comparison
#   7. apeglm log2 fold-change shrinkage
#   8. PCA plot + DE (volcano) plot
#   9. Export complete shrunken results table
#  10. Save fitted object + sessionInfo()
#  11. AI use is documented in week5_AI_verification_log.md
#  12. Interpretation written in week5_interpretation.md
#
# Run from any directory:
#   Rscript week5_deseq2_analysis.R
# Input files are expected one level up in the "for_student" folder.
# ================================================================

# ---- 0. Working directory and paths ---------------------------------
args <- commandArgs(trailingOnly = FALSE)
file_arg <- args[grep("^--file=", args)]
if (length(file_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
} else {
  script_dir <- getwd()
}
out_dir  <- script_dir                                  # deliverables land here
input_dir <- normalizePath(file.path(script_dir, "..", "for_student"))

count_file     <- file.path(input_dir, "Week5_Homework_Count_Matrix.csv")
metadata_file  <- file.path(input_dir, "Week5_Homework_Sample_Metadata.csv")

# ---- 1. Load packages -------------------------------------------------
suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggrepel)
})

# ---- 2. Import ---------------------------------------------------------
counts <- read.csv(count_file, row.names = 1, check.names = FALSE)
coldata <- read.csv(metadata_file, row.names = 1, check.names = FALSE)

cat("Count matrix dimensions:", nrow(counts), "genes x", ncol(counts), "samples\n")
cat("Metadata dimensions:", nrow(coldata), "samples x", ncol(coldata), "columns\n")

# ---- 3. Mandatory validation ------------------------------------------
# 3.1 Non-negative integers
stopifnot("counts must be numeric" = is.numeric(as.matrix(counts)))
stopifnot("counts must be non-negative" = all(counts >= 0))
stopifnot("counts must be integers" = all(as.matrix(counts) == round(as.matrix(counts))))

# 3.2 Columns of count matrix exactly match metadata row names (same order)
stopifnot("ncol(counts) != nrow(coldata)" = ncol(counts) == nrow(coldata))
stopifnot("column/row name mismatch" = identical(colnames(counts), rownames(coldata)))

# 3.3 No duplicated sample IDs
stopifnot("duplicated sample ids" = !anyDuplicated(rownames(coldata)))

# 3.4 Factors with expected levels
coldata$condition <- relevel(factor(coldata$condition), ref = "control")
coldata$batch     <- factor(coldata$batch)
stopifnot("unexpected condition levels" =
            setequal(levels(coldata$condition), c("control", "treated")))
stopifnot("unexpected batch levels" =
            setequal(levels(coldata$batch), c("A", "B", "C")))

# 3.5 control is the reference level
stopifnot("control must be reference" = levels(coldata$condition)[1] == "control")

cat("\nDesign table (batch x condition):\n")
print(table(coldata$batch, coldata$condition))
cat("\nLibrary sizes (total counts per sample):\n")
print(sort(colSums(counts)))

# ---- 4. Construct DESeq2 object ---------------------------------------
# Batch is included in the design because the 12 samples come from three
# balanced batches (A/B/C); DESeq2 then models batch effects and estimates
# the treatment effect conditional on batch, avoiding batch-induced
# inflation of the condition contrast.
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = coldata,
  design    = ~ batch + condition
)

# ---- 5. Pre-filter ------------------------------------------------------
# Retain genes with at least 10 counts in at least 3 samples.
keep <- rowSums(counts(dds) >= 10) >= 3
cat("\nGenes before filtering:", nrow(dds), "\n")
dds <- dds[keep, ]
cat("Genes after filtering:", nrow(dds), "\n")

# ---- 6. Fit model --------------------------------------------------------
dds <- DESeq(dds)

coef_names <- resultsNames(dds)
cat("\nresultsNames(dds):\n")
print(coef_names)

# Confirm the exact treated-vs-control coefficient name
target_coef <- "condition_treated_vs_control"
if (!target_coef %in% coef_names) {
  stop("Expected coefficient was not found. Inspect resultsNames(dds) and update target_coef.")
}

# ---- 7. Extract contrast and shrink ------------------------------------
# Alpha controls the FDR cut-off used by DESeq2's independent filtering;
# final significance uses padj < 0.05 and |log2FC| >= 1.
res <- results(dds, contrast = c("condition", "treated", "control"), alpha = 0.05)

res_shrunk <- lfcShrink(dds, coef = target_coef, type = "apeglm")

# ---- 8. Assemble results table -------------------------------------------
res_df <- as.data.frame(res_shrunk) |>
  rownames_to_column("gene_id") |>
  mutate(
    significant = !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1,
    direction = case_when(
      significant & log2FoldChange > 0 ~ "Up in treated",
      significant & log2FoldChange < 0 ~ "Down in treated",
      TRUE ~ "Not significant"
    )
  ) |>
  arrange(padj)

write.csv(res_df, file.path(out_dir, "week5_deseq2_results.csv"), row.names = FALSE)

cat("\nSignificant genes (padj < 0.05, |log2FC| >= 1):", sum(res_df$significant), "\n")
print(table(res_df$direction))

# ---- 9. PCA --------------------------------------------------------------
# After pre-filtering, 989 genes remain (< nsub = 1000 required by the fast
# vst() approximation), so we call varianceStabilizingTransformation directly,
# which computes the same VST without the subset approximation.
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

pca_df <- plotPCA(vsd, intgroup = c("condition", "batch"), returnData = TRUE)
percent_var <- round(100 * attr(pca_df, "percentVar"))

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2,
                            color = condition, shape = batch, label = name)) +
  geom_point(size = 4) +
  geom_text_repel(size = 3, max.overlaps = Inf) +
  labs(
    title = "Week 5 RNA-seq PCA (VST)",
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance")
  ) +
  theme_bw(base_size = 12)

ggsave(file.path(out_dir, "week5_pca.png"), p_pca, width = 7, height = 5, dpi = 300)

# ---- 10. DE plot (volcano) ------------------------------------------------
plot_df <- res_df |>
  mutate(neg_log10_padj = -log10(pmax(padj, 1e-300)))

p_volcano <- ggplot(plot_df, aes(x = log2FoldChange, y = neg_log10_padj, color = direction)) +
  geom_point(alpha = 0.7, size = 1.8) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_color_manual(values = c(
    "Up in treated" = "#C0392B",
    "Down in treated" = "#2F6DB3",
    "Not significant" = "grey70"
  )) +
  labs(
    title = "Differential expression: treated versus control",
    x = "Shrunken log2 fold change (apeglm)",
    y = "-log10 adjusted p value",
    color = NULL
  ) +
  theme_bw(base_size = 12)

ggsave(file.path(out_dir, "week5_de_plot.png"), p_volcano, width = 7, height = 5, dpi = 300)

# ---- 11. Save reproducibility files ----------------------------------------
saveRDS(dds, file.path(out_dir, "week5_deseq2_object.rds"))

sink(file.path(out_dir, "session_info.txt"))
sessionInfo()
sink()

# ---- 12. Summary printed for interpretation --------------------------------
cat("\n=== Summary for interpretation ===\n")
cat("Comparison: treated versus control (design ~ batch + condition)\n")
cat("Significant up in treated:", sum(res_df$direction == "Up in treated"), "\n")
cat("Significant down in treated:", sum(res_df$direction == "Down in treated"), "\n")
cat("Top 10 by padj:\n")
print(head(res_df, 10))
cat("\nAnalysis completed. Deliverables written to:", out_dir, "\n")
