# ============================================================
# MDD Gene Expression Analysis — GSE98793
# 03: Differential Expression Analysis with limma
# ============================================================

library(limma)
library(ggplot2)

# ---- 1. Load preprocessed data ----
load("data/processed/MDD_batch_corrected.RData")

cat("Expression matrix dimensions:", dim(expr_filtered), "\n")
cat("Group counts:\n")
print(table(pheno_factor))

# ---- 2. Create design matrix ----
batch <- pheno_data$`batch:ch1`

design <- model.matrix(~ pheno_factor + batch)


# ---- 3. Fit limma model ----
#design matrix constructed in preprocessing step:~pheno_factor
fit <- lmFit(expr_corrected, design)
fit <- eBayes(fit)

# ---- 4. Extract all differential expression results ----
de_results <- topTable(
  fit,
  coef = "pheno_factorMDD",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

# ---- 5. Add feature annotations ----
de_results$ID <- rownames(de_results)

annotation_cols <- c("ID", "Gene Symbol", "Gene Title", "ENTREZ_GENE_ID")
annotation_df <- feature_data[, annotation_cols]

de_results_annotated <- merge(
  de_results,
  annotation_df,
  by = "ID",
  all.x = TRUE,
  sort = FALSE
)

# ---- 6. Reorder columns using base R ----
de_results_annotated <- de_results_annotated[, c(
  "ID",
  "Gene Symbol",
  "Gene Title",
  "ENTREZ_GENE_ID",
  "logFC",
  "AveExpr",
  "t",
  "P.Value",
  "adj.P.Val",
  "B"
)]
# ---- 7. Filter significant genes ----
sig_genes <- de_results_annotated[
  !is.na(de_results_annotated$`Gene Symbol`) &
    de_results_annotated$adj.P.Val < 0.05 &
    abs(de_results_annotated$logFC) > 0.3,
]
cat("\nNumber of significant genes:", nrow(sig_genes), "\n")

# ---- 8. Create volcano plot labels ----
de_results_annotated$significance <- "Not Significant"
de_results_annotated$significance[
  de_results_annotated$adj.P.Val < 0.05 & de_results_annotated$logFC > 0.3
] <- "Up in MDD"
de_results_annotated$significance[
  de_results_annotated$adj.P.Val < 0.05 & de_results_annotated$logFC < -0.3
] <- "Down in MDD"

# ---- 9. Volcano plot ----
volcano_plot <- ggplot(de_results_annotated, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = significance), alpha = 0.6, size = 1.2) +
  scale_color_manual(values = c(
    "Up in MDD" = "red",
    "Down in MDD" = "blue",
    "Not Significant" = "grey"
  )) +
  geom_vline(xintercept = c(-0.3, 0.3), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "Volcano Plot: MDD vs HC",
    x = "log2 Fold Change",
    y = "-log10 Adjusted P-value"
  )

# ---- 10. Save results ----
dir.create("results", recursive = TRUE, showWarnings = FALSE)

write.csv(de_results_annotated,
          "results/MDD_DE_results_full.csv",
          row.names = FALSE)

write.csv(sig_genes,
          "results/MDD_DE_genes_sig.csv",
          row.names = FALSE)

ggsave("results/MDD_volcano_plot.png",
       volcano_plot)

# ---- 11. Print top results ----
cat("\nTop 10 DE genes:\n")
print(head(de_results_annotated, 10))
