# ============================================================
# MDD Gene Expression Analysis — GSE98793
# 02: Differential Expression Analysis with limma
# ============================================================

library(limma)
library(ggplot2)

# ---- 1. Load preprocessed data ----
load("~/Desktop/MDD_preprocessed.RData")

cat("Expression matrix dimensions:", dim(expr_filtered), "\n")
cat("Group counts:\n")
print(table(pheno_factor))

# ---- 2. Fit limma model ----
fit <- lmFit(expr_filtered, design)
fit <- eBayes(fit)

# ---- 3. Extract all differential expression results ----
de_results <- topTable(
  fit,
  coef = "pheno_factorMDD",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

# ---- 4. Add feature annotations ----
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

# ---- 5. Reorder columns using base R ----
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

# ---- 6. Save full results ----
write.csv(
  de_results_annotated,
  file = "~/Desktop/MDD_DE_results_full.csv",
  row.names = FALSE
)

cat("\nSaved full DE results: ~/Desktop/MDD_DE_results_full.csv\n")

# ---- 7. Filter significant genes ----
sig_genes <- de_results_annotated[
  !is.na(de_results_annotated$`Gene Symbol`) &
    de_results_annotated$adj.P.Val < 0.05 &
    abs(de_results_annotated$logFC) > 1,
]

write.csv(
  sig_genes,
  file = "~/Desktop/MDD_DE_genes_sig.csv",
  row.names = FALSE
)

cat("Saved significant DE genes: ~/Desktop/MDD_DE_genes_sig.csv\n")
cat("\nNumber of significant genes:", nrow(sig_genes), "\n")

# ---- 8. Create volcano plot labels ----
de_results_annotated$significance <- "Not Significant"
de_results_annotated$significance[
  de_results_annotated$adj.P.Val < 0.05 & de_results_annotated$logFC > 1
] <- "Up in MDD"
de_results_annotated$significance[
  de_results_annotated$adj.P.Val < 0.05 & de_results_annotated$logFC < -1
] <- "Down in MDD"

# ---- 9. Volcano plot ----
volcano_plot <- ggplot(de_results_annotated, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = significance), alpha = 0.6, size = 1.2) +
  scale_color_manual(values = c(
    "Up in MDD" = "red",
    "Down in MDD" = "blue",
    "Not Significant" = "grey"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "Volcano Plot: MDD vs HC",
    x = "log2 Fold Change",
    y = "-log10 Adjusted P-value"
  )

ggsave(
  filename = "~/Desktop/MDD_volcano_plot.png",
  plot = volcano_plot,
  width = 8,
  height = 6,
  dpi = 300
)

cat("Saved volcano plot: ~/Desktop/MDD_volcano_plot.png\n")

# ---- 10. Print top results ----
cat("\nTop 10 DE genes:\n")
print(head(de_results_annotated, 10))
