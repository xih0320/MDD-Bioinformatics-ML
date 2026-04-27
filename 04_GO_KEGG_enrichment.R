# ===============================
# 04_GO_KEGG_enrichment.R
# GO and KEGG Enrichment Analysis
# ===============================

# 1. Load required packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)

dir.create("results", showWarnings = FALSE)
dir.create("plots", showWarnings = FALSE)

# 2. Define file paths
sig_gene_file <- "results/MDD_DE_genes_sig.csv"

go_result_file <- "results/GO_enrichment_results.csv"
kegg_result_file <- "results/KEGG_enrichment_results.csv"

mapped_gene_file <- "results/mapped_gene_symbols_to_entrez.csv"

go_dotplot_file <- "plots/GO_dotplot.png"
kegg_dotplot_file <- "plots/KEGG_dotplot.png"

# 3. Load significant DE genes
de_genes <- read.csv(sig_gene_file, check.names = FALSE)

# 4. Check column names
print(colnames(de_genes))
head(de_genes)

# 5. Extract gene symbols
gene_symbols <- de_genes$`Gene Symbol`

# Remove missing or empty gene symbols
gene_symbols <- gene_symbols[!is.na(gene_symbols)]
gene_symbols <- gene_symbols[gene_symbols != ""]
gene_symbols <- unique(gene_symbols)

cat("Number of unique gene symbols:", length(gene_symbols), "\n")

# 6. Convert gene symbols to Entrez IDs
gene_df <- bitr(
  gene_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

gene_ids <- unique(gene_df$ENTREZID)

cat("Number of successfully mapped Entrez IDs:", length(gene_ids), "\n")

# Save mapped gene ID table
write.csv(
  gene_df,
  mapped_gene_file,
  row.names = FALSE
)

# ===============================
# 7. GO Enrichment Analysis
# ===============================

go_result <- enrichGO(
  gene          = gene_ids,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

go_result_df <- as.data.frame(go_result)

write.csv(
  go_result_df,
  go_result_file,
  row.names = FALSE
)

cat("GO enrichment results saved to:", go_result_file, "\n")
cat("Number of significant GO terms:", nrow(go_result_df), "\n")

if (nrow(go_result_df) > 0) {
  png(go_dotplot_file, width = 1000, height = 800)
  print(dotplot(go_result, showCategory = 10))
  dev.off()
  cat("GO dotplot saved to:", go_dotplot_file, "\n")
} else {
  cat("No significant GO terms found. GO dotplot was not generated.\n")
}

# ===============================
# 8. KEGG Enrichment Analysis
# ===============================

kegg_result <- enrichKEGG(
  gene          = gene_ids,
  organism      = "hsa",
  keyType       = "kegg",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH"
)

kegg_result <- setReadable(
  kegg_result,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
)

kegg_result_df <- as.data.frame(kegg_result)

write.csv(
  kegg_result_df,
  kegg_result_file,
  row.names = FALSE
)

cat("KEGG enrichment results saved to:", kegg_result_file, "\n")
cat("Number of significant KEGG pathways:", nrow(kegg_result_df), "\n")

if (nrow(kegg_result_df) > 0) {
  png(kegg_dotplot_file, width = 1000, height = 800)
  print(dotplot(kegg_result, showCategory = 10))
  dev.off()
  cat("KEGG dotplot saved to:", kegg_dotplot_file, "\n")
} else {
  cat("No significant KEGG pathways found. KEGG dotplot was not generated.\n")
}

# ===============================
# 9. Print top results
# ===============================

cat("\nTop GO enrichment results:\n")
print(head(go_result_df, 10))

cat("\nTop KEGG enrichment results:\n")
print(head(kegg_result_df, 10))

cat("\nGO and KEGG enrichment analysis completed successfully.\n")