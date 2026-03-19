# ============================================================
# MDD Gene Expression Analysis
# 02: Differential Expression with t-test
# ============================================================
# Requires: GxS_covfilter, pheno_factor from 01_data_preprocessing.R

library(ggplot2)
library(dplyr)

# ---- t-test for a single gene ----
myrow = 2
mygene = rownames(GxS_covfilter)[myrow]
cat("Testing gene:", mygene, "\n")

# Formula interface (cleaner for large-scale use)
t_result = t.test(GxS_covfilter[myrow, ] ~ pheno_factor)
cat("p-value:", t_result$p.value, "\n")

# Boxplot for single gene
mygene_data_df = data.frame(
  gene = GxS_covfilter[myrow, ],
  phenotype = pheno_factor
)

p = ggplot(mygene_data_df, aes(x = phenotype, y = gene, fill = phenotype)) +
  stat_boxplot(geom = 'errorbar') +
  geom_boxplot() +
  xlab("MDD vs HC") +
  ylab(mygene) +
  theme_minimal()
print(p)

# ---- t-test across all genes ----
ttest_fn = function(i) {
  mygene  = rownames(GxS_covfilter)[i]
  t_result = t.test(GxS_covfilter[i, ] ~ pheno_factor)
  tstat    = t_result$statistic
  pval     = t_result$p.value
  c(mygene, tstat, pval)
}

# Apply to all genes using sapply
ttest_allgene_mat = t(sapply(1:nrow(GxS_covfilter), ttest_fn))
ttest_allgene_df  = data.frame(ttest_allgene_mat)
colnames(ttest_allgene_df) = c("gene", "tstat", "pval")

# Sort by p-value
ttest_allgene_sorted = ttest_allgene_df %>%
  mutate_at("pval", as.character) %>%
  mutate_at("pval", as.numeric) %>%
  arrange(pval)

cat("\nTop 10 differentially expressed genes:\n")
print(ttest_allgene_sorted[1:10, ])

# ---- Boxplot for top gene ----
top_gene = as.character(ttest_allgene_sorted[1, 1])
myrow    = which(rownames(GxS_covfilter) == top_gene)

top_gene_df = data.frame(
  gene = GxS_covfilter[myrow, ],
  phenotype = pheno_factor
)

p2 = ggplot(top_gene_df, aes(x = phenotype, y = gene, fill = phenotype)) +
  stat_boxplot(geom = 'errorbar') +
  geom_boxplot() +
  xlab("MDD vs HC") +
  ylab(top_gene) +
  ggtitle(paste("Top DE gene:", top_gene)) +
  theme_minimal()
print(p2)

# ---- Export top 200 genes for pathway enrichment (MSigDB Reactome) ----
top_cutoff = 200
top_genes  = as.character(ttest_allgene_sorted[1:top_cutoff, 1])
write.table(top_genes, sep = "\t", file = "top200_DE_genes.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)
cat("\nTop 200 genes written to top200_DE_genes.txt\n")
