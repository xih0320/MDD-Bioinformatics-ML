# ===============================
# 04_GO_KEGG_GSEA_enrichment.R
# GO, KEGG, and GSEA-KEGG Enrichment Analysis
# ===============================

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)

set.seed(123)

dir.create("results", showWarnings = FALSE)
dir.create("plots", showWarnings = FALSE)

# ===============================
# 1. Define file paths
# ===============================

sig_gene_file  <- "results/MDD_DE_genes_sig.csv"
full_de_file   <- "results/MDD_DE_results_full.csv"

go_result_file        <- "results/GO_enrichment_results.csv"
kegg_result_file      <- "results/KEGG_enrichment_results.csv"
gsea_kegg_result_file <- "results/GSEA_KEGG_results.csv"

mapped_sig_gene_file  <- "results/mapped_sig_genes_to_entrez.csv"
mapped_full_gene_file <- "results/mapped_full_genes_to_entrez.csv"

go_dotplot_file        <- "plots/GO_dotplot.png"
kegg_dotplot_file      <- "plots/KEGG_dotplot.png"
gsea_kegg_dotplot_file <- "plots/GSEA_KEGG_dotplot.png"

# ===============================
# 2. Load data
# ===============================

sig_de_genes <- read.csv(sig_gene_file, check.names = FALSE)
full_de <- read.csv(full_de_file, check.names = FALSE)

cat("Significant DE gene file columns:\n")
print(colnames(sig_de_genes))

cat("Full DE result file columns:\n")
print(colnames(full_de))

# ===============================
# 3. Prepare significant genes for GO / KEGG ORA
# ===============================

sig_gene_symbols <- sig_de_genes$`Gene Symbol`

sig_gene_symbols <- sig_gene_symbols[!is.na(sig_gene_symbols)]
sig_gene_symbols <- sig_gene_symbols[sig_gene_symbols != ""]
sig_gene_symbols <- unique(sig_gene_symbols)

cat("Number of unique significant gene symbols:", length(sig_gene_symbols), "\n")

# Convert significant gene symbols to Entrez IDs
sig_gene_df <- bitr(
  sig_gene_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

sig_gene_ids <- unique(sig_gene_df$ENTREZID)

cat("Number of successfully mapped significant Entrez IDs:", length(sig_gene_ids), "\n")

write.csv(
  sig_gene_df,
  mapped_sig_gene_file,
  row.names = FALSE
)

# ===============================
# 4. GO Enrichment Analysis
# ===============================

go_result <- enrichGO(
  gene          = sig_gene_ids,
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
# 5. Traditional KEGG Enrichment Analysis
# ===============================

kegg_result <- enrichKEGG(
  gene          = sig_gene_ids,
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
# 6. Prepare full ranked gene list for GSEA-KEGG
# ===============================

full_de_clean <- full_de %>%
  filter(
    !is.na(`Gene Symbol`),
    `Gene Symbol` != "",
    !is.na(logFC)
  )

full_gene_df <- bitr(
  full_de_clean$`Gene Symbol`,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

write.csv(
  full_gene_df,
  mapped_full_gene_file,
  row.names = FALSE
)

full_de_mapped <- full_de_clean %>%
  inner_join(
    full_gene_df,
    by = c("Gene Symbol" = "SYMBOL")
  )

# If multiple probes map to the same gene,
# keep the probe with the largest absolute logFC
full_de_mapped <- full_de_mapped %>%
  group_by(ENTREZID) %>%
  slice_max(order_by = abs(logFC), n = 1, with_ties = FALSE) %>%
  ungroup()

gene_list <- full_de_mapped$logFC
names(gene_list) <- full_de_mapped$ENTREZID

gene_list <- sort(gene_list, decreasing = TRUE)

cat("Number of genes used for GSEA-KEGG:", length(gene_list), "\n")

# ===============================
# 7. GSEA-KEGG Analysis
# ===============================

gsea_kegg_result <- gseKEGG(
  geneList     = gene_list,
  organism     = "hsa",
  keyType      = "ncbi-geneid",
  pvalueCutoff = 0.05,
  verbose      = FALSE
)

gsea_kegg_result <- setReadable(
  gsea_kegg_result,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
)

gsea_kegg_df <- as.data.frame(gsea_kegg_result)

write.csv(
  gsea_kegg_df,
  gsea_kegg_result_file,
  row.names = FALSE
)

cat("GSEA-KEGG results saved to:", gsea_kegg_result_file, "\n")
cat("Number of significant GSEA-KEGG pathways:", nrow(gsea_kegg_df), "\n")

if (nrow(gsea_kegg_df) > 0) {
  png(gsea_kegg_dotplot_file, width = 1000, height = 800)
  print(dotplot(gsea_kegg_result, showCategory = 10))
  dev.off()
  cat("GSEA-KEGG dotplot saved to:", gsea_kegg_dotplot_file, "\n")
} else {
  cat("No significant GSEA-KEGG pathways found.\n")
}

# ===============================
# 8. Print top results
# ===============================

cat("\nTop GO enrichment results:\n")
print(head(go_result_df, 10))

cat("\nTop traditional KEGG enrichment results:\n")
print(head(kegg_result_df, 10))

cat("\nTop GSEA-KEGG results:\n")
print(head(gsea_kegg_df, 10))

cat("\nGO, KEGG, and GSEA-KEGG enrichment analysis completed successfully.\n")