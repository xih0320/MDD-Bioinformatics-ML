# ============================================================
# MDD Gene Expression Analysis
# 03: Feature Selection — Logistic Regression & LASSO
# ============================================================
# Requires: GxS_covfilter, pheno_factor from 01_data_preprocessing.R

library(dplyr)
library(ggplot2)
library(glmnet)

# Recode phenotype: HC = 0 (reference), MDD = 1
pheno_factor_relevel = relevel(pheno_factor, "HC")
levels(pheno_factor_relevel)[levels(pheno_factor_relevel) == "MDD"] = 1
levels(pheno_factor_relevel)[levels(pheno_factor_relevel) == "HC"]  = 0

# ---- Univariate logistic regression for all genes ----
lr.fn = function(i) {
  gene      = rownames(GxS_covfilter)[i]
  gene_expr = GxS_covfilter[i, ]
  gene_fit  = glm(pheno_factor_relevel ~ gene_expr, family = binomial)
  coeff_mat = coef(summary(gene_fit))
  b1        = coeff_mat[2, 1]
  b1_pval   = coeff_mat[2, 4]
  c(gene, b1, b1_pval)
}

num.genes      = nrow(GxS_covfilter)
lr_results_mat = matrix(0, nrow = num.genes, ncol = 3)

for (i in 1:num.genes) {
  lr_results_mat[i, ] = lr.fn(i)
}

lr_results_df = data.frame(lr_results_mat)
colnames(lr_results_df) = c("gene", "b1", "pval")

lr_results_sorted = lr_results_df %>%
  mutate_at("pval", as.character) %>%
  mutate_at("pval", as.numeric) %>%
  arrange(pval)

cat("Top 10 genes by logistic regression p-value:\n")
print(lr_results_sorted[1:10, ])

# FDR-adjusted p-values
lr_p_adj = p.adjust(lr_results_sorted$pval, method = "BH")
lr_results_sorted$p_adj = lr_p_adj
cat("\nTop genes with FDR < 0.05:\n")
print(lr_results_sorted[lr_results_sorted$p_adj < 0.05, ])

# ---- LASSO with cross-validation (glmnet) ----
# Transpose so subjects are rows
X = t(GxS_covfilter)
y = as.numeric(pheno_factor_relevel) - 1  # 0/1

set.seed(42)
glmnet_model = cv.glmnet(
  x            = X,
  y            = y,
  alpha        = 1,          # alpha=1 -> LASSO
  family       = "binomial",
  type.measure = "class",
  keep         = TRUE
)

plot(glmnet_model)
cat("\nOptimal lambda (min CV error):", glmnet_model$lambda.min, "\n")

# ---- FIX: extract non-zero LASSO coefficients ----
glmnet_coeffs = predict(glmnet_model, type = "coefficients",
                        s = glmnet_model$lambda.min)
coeff_df = as.data.frame(as.matrix(glmnet_coeffs))
coeff_df$features = rownames(coeff_df)
colnames(coeff_df)[1] = "coef"
top_glmnet = coeff_df[coeff_df$coef != 0 &
                      coeff_df$features != "(Intercept)", ]

cat("\nGenes selected by LASSO:", nrow(top_glmnet), "\n")
print(top_glmnet)

# Classification accuracy (in-sample)
glmnet_predicted = predict(glmnet_model,
                           s    = glmnet_model$lambda.min,
                           type = "class",
                           newx = X)

glmnet_accuracy = mean(factor(glmnet_predicted) == pheno_factor_relevel)
cat("\nLASSO Classification Accuracy:", round(glmnet_accuracy * 100, 2), "%\n")

# Confusion matrix
cat("\nConfusion matrix:\n")
print(table(Predicted = glmnet_predicted, Actual = y))
