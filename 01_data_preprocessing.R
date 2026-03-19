# ============================================================
# MDD Gene Expression Analysis
# 01: Data Loading, Normalization, and Preprocessing
# ============================================================

# Load gene expression data (RNA-seq counts per million)
load("sense.filtered.cpm.Rdata")

# Load clinical/demographic data
subject.attrs = read.csv("Demographic_symptom.csv", stringsAsFactors = FALSE)

library(dplyr)

# Match subject IDs between expression and phenotype data
phenos_df = subject.attrs %>%
  filter(X %in% colnames(sense.filtered.cpm)) %>%
  dplyr::select(X, Diag)

mddPheno = as.factor(phenos_df$Diag)

# Quantile normalization + log2 transformation
library(preprocessCore)

mddExprData_quantile = normalize.quantiles(sense.filtered.cpm)
mddExprData_quantileLog2 = log2(mddExprData_quantile)

# Attach gene and phenotype names
colnames(mddExprData_quantileLog2) = mddPheno
rownames(mddExprData_quantileLog2) = rownames(sense.filtered.cpm)

# Total genes before filtering
cat("Genes before filtering:", nrow(mddExprData_quantileLog2), "\n")

# ---- Coefficient of Variation (CoV) Filter ----
# CoV = sd(x) / abs(mean(x))
# Low CoV = experimental effect large relative to measurement noise
CoV_values = apply(mddExprData_quantileLog2, 1, function(x) { sd(x) / abs(mean(x)) })

# Remove zero-variance genes
sd_values = apply(mddExprData_quantileLog2, 1, function(x) { sd(x) })
cat("Zero-variance gene:", rownames(mddExprData_quantileLog2)[sd_values == 0], "\n")

# Apply filter
GxS_covfilter = mddExprData_quantileLog2[CoV_values < .045 & sd_values > 0, ]
cat("Genes after CoV filter:", nrow(GxS_covfilter), "\n")
cat("Subjects:", ncol(GxS_covfilter), "\n")

# Phenotype factor for downstream analysis
pheno_factor = as.factor(colnames(GxS_covfilter))
cat("Phenotype levels:", levels(pheno_factor), "\n")
