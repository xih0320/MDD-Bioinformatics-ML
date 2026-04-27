# ==============================
# 02: QC & PCA Analysis
# ==============================

library(ggplot2)
library(sva)

# ---- 1. Load processed data ----
load("data/processed/MDD_preprocessed.RData")

table(pheno_factor)
colnames(pheno_data)

# ---- 2. Define batch variable ----
batch <- pheno_data$`batch:ch1`

# ---- 3. PCA before batch correction ----
expr_t_before <- t(expr_filtered)

pca_before <- prcomp(expr_t_before, scale. = TRUE)

pca_before_df <- data.frame(
  PC1 = pca_before$x[, 1],
  PC2 = pca_before$x[, 2],
  Group = pheno_factor,
  Batch = batch
)

p_before_group <- ggplot(pca_before_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3, alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "PCA Before Batch Correction",
    subtitle = "Colored by phenotype (HC vs MDD)",
    x = "PC1",
    y = "PC2"
  )

p_before_batch <- ggplot(pca_before_df, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size = 3, alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "PCA Before Batch Correction",
    subtitle = "Colored by batch",
    x = "PC1",
    y = "PC2"
  )

print(p_before_group)
print(p_before_batch)

# ---- 4. Batch correction with ComBat ----
expr_corrected <- ComBat(
  dat = expr_filtered,
  batch = batch,
  par.prior = TRUE,
  prior.plots = FALSE
)

# ---- 5. PCA after batch correction ----
expr_t_after <- t(expr_corrected)

pca_after <- prcomp(expr_t_after, scale. = TRUE)

pca_after_df <- data.frame(
  PC1 = pca_after$x[, 1],
  PC2 = pca_after$x[, 2],
  Group = pheno_factor,
  Batch = batch
)

p_after_group <- ggplot(pca_after_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3, alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "PCA After Batch Correction",
    subtitle = "Colored by phenotype (HC vs MDD)",
    x = "PC1",
    y = "PC2"
  )

p_after_batch <- ggplot(pca_after_df, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size = 3, alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "PCA After Batch Correction",
    subtitle = "Colored by batch",
    x = "PC1",
    y = "PC2"
  )

print(p_after_group)
print(p_after_batch)

# ---- 6. Save plots ----
dir.create("results/plots", recursive = TRUE, showWarnings = FALSE)

ggsave("results/plots/PCA_before_group.png", p_before_group, width = 6, height = 5, dpi = 300)
ggsave("results/plots/PCA_before_batch.png", p_before_batch, width = 6, height = 5, dpi = 300)
ggsave("results/plots/PCA_after_group.png", p_after_group, width = 6, height = 5, dpi = 300)
ggsave("results/plots/PCA_after_batch.png", p_after_batch, width = 6, height = 5, dpi = 300)

# ---- 7. Save corrected expression data ----
save(
  expr_corrected,
  pheno_factor,
  pheno_data,
  feature_data,
  design,
  file = "data/processed/MDD_batch_corrected.RData"
)

# ---- 8. Check variance explained ----
variance_before <- summary(pca_before)$importance[2, 1:2]
variance_after <- summary(pca_after)$importance[2, 1:2]

cat("Variance explained before batch correction:\n")
print(variance_before)

cat("Variance explained after batch correction:\n")
print(variance_after)