# ============================================================
# MDD Gene Expression Analysis
# 06: Machine Learning Model Comparison
# LASSO vs Random Forest vs XGBoost
# ============================================================

library(glmnet)
library(ranger)
library(xgboost)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)

# ---- 1. Create output folders ----
dir.create("results", showWarnings = FALSE)
dir.create("plots", showWarnings = FALSE)

# ---- 2. Load batch-corrected expression data and phenotype labels ----
expr_matrix <- readRDS("results/batch_corrected_expression_matrix.rds")
pheno_factor <- readRDS("results/pheno_factor.rds")

cat("Expression matrix dimensions:\n")
print(dim(expr_matrix))

cat("Group distribution:\n")
print(table(pheno_factor))

# ---- 3. Load differential expression results ----
de_results <- read.csv("results/MDD_DE_results_full.csv", check.names = FALSE)

# Use top 200 ranked probes based on adjusted p-value
top_features <- de_results %>%
  filter(!is.na(ID)) %>%
  arrange(adj.P.Val) %>%
  pull(ID) %>%
  unique() %>%
  head(200)

common_features <- intersect(top_features, rownames(expr_matrix))

cat("Number of top DE features used for ML:", length(common_features), "\n")

# Save selected feature list
selected_feature_table <- de_results %>%
  filter(ID %in% common_features) %>%
  select(ID, `Gene Symbol`, `Gene Title`, ENTREZ_GENE_ID, logFC, P.Value, adj.P.Val)

write.csv(
  selected_feature_table,
  "results/ML_top200_selected_features.csv",
  row.names = FALSE
)

# ---- 4. Prepare ML matrix ----
X <- t(expr_matrix[common_features, ])
y <- pheno_factor

# Make sure labels are clean
y <- factor(y, levels = c("HC", "MDD"))

cat("ML input dimensions:\n")
print(dim(X))
print(table(y))

# ---- 5. Train/test split ----
set.seed(42)

train_index <- createDataPartition(y, p = 0.8, list = FALSE)

X_train <- X[train_index, ]
X_test  <- X[-train_index, ]

y_train <- y[train_index]
y_test  <- y[-train_index]

cat("Training set:\n")
print(table(y_train))

cat("Test set:\n")
print(table(y_test))

# Numeric labels for glmnet and xgboost
y_train_num <- ifelse(y_train == "MDD", 1, 0)
y_test_num  <- ifelse(y_test == "MDD", 1, 0)

# ============================================================
# 6. LASSO Logistic Regression
# ============================================================

set.seed(42)

lasso_cv <- cv.glmnet(
  x = X_train,
  y = y_train_num,
  family = "binomial",
  alpha = 1,
  type.measure = "auc",
  nfolds = 5
)

lasso_prob <- as.numeric(
  predict(lasso_cv, newx = X_test, s = "lambda.min", type = "response")
)

lasso_pred <- factor(
  ifelse(lasso_prob >= 0.5, "MDD", "HC"),
  levels = c("HC", "MDD")
)

lasso_cm <- confusionMatrix(lasso_pred, y_test, positive = "MDD")
lasso_roc <- roc(y_test_num, lasso_prob, quiet = TRUE)
lasso_auc <- as.numeric(auc(lasso_roc))

cat("\nLASSO Results:\n")
print(lasso_cm)
cat("LASSO AUC:", round(lasso_auc, 4), "\n")

# Save non-zero LASSO coefficients
lasso_coef <- coef(lasso_cv, s = "lambda.min")
lasso_coef_df <- data.frame(
  Feature = rownames(lasso_coef),
  Coefficient = as.numeric(lasso_coef)
)

lasso_selected <- lasso_coef_df %>%
  filter(Coefficient != 0, Feature != "(Intercept)") %>%
  left_join(
    selected_feature_table,
    by = c("Feature" = "ID")
  ) %>%
  arrange(desc(abs(Coefficient)))

write.csv(
  lasso_selected,
  "results/LASSO_selected_features.csv",
  row.names = FALSE
)

# ============================================================
# 7. Random Forest
# ============================================================

rf_train_df <- data.frame(X_train)
rf_train_df$Group <- y_train

rf_test_df <- data.frame(X_test)

set.seed(42)

rf_model <- ranger(
  Group ~ .,
  data = rf_train_df,
  num.trees = 1000,
  probability = TRUE,
  importance = "permutation",
  seed = 42
)

rf_prob <- predict(rf_model, data = rf_test_df)$predictions[, "MDD"]

rf_pred <- factor(
  ifelse(rf_prob >= 0.5, "MDD", "HC"),
  levels = c("HC", "MDD")
)

rf_cm <- confusionMatrix(rf_pred, y_test, positive = "MDD")
rf_roc <- roc(y_test_num, rf_prob, quiet = TRUE)
rf_auc <- as.numeric(auc(rf_roc))

cat("\nRandom Forest Results:\n")
print(rf_cm)
cat("Random Forest AUC:", round(rf_auc, 4), "\n")

# Save RF feature importance
rf_importance <- data.frame(
  Feature = names(rf_model$variable.importance),
  Importance = as.numeric(rf_model$variable.importance)
) %>%
  arrange(desc(Importance)) %>%
  left_join(
    selected_feature_table,
    by = c("Feature" = "ID")
  )

write.csv(
  rf_importance,
  "results/RandomForest_feature_importance.csv",
  row.names = FALSE
)

# ============================================================
# 8. XGBoost
# ============================================================

dtrain <- xgb.DMatrix(data = X_train, label = y_train_num)
dtest  <- xgb.DMatrix(data = X_test, label = y_test_num)

set.seed(42)

xgb_model <- xgb.train(
  params = list(
    objective = "binary:logistic",
    eval_metric = "auc",
    eta = 0.05,
    max_depth = 3,
    subsample = 0.8,
    colsample_bytree = 0.8
  ),
  data = dtrain,
  nrounds = 100,
  verbose = 0
)

xgb_prob <- predict(xgb_model, dtest)

xgb_pred <- factor(
  ifelse(xgb_prob >= 0.5, "MDD", "HC"),
  levels = c("HC", "MDD")
)

xgb_cm <- confusionMatrix(xgb_pred, y_test, positive = "MDD")
xgb_roc <- roc(y_test_num, xgb_prob, quiet = TRUE)
xgb_auc <- as.numeric(auc(xgb_roc))

cat("\nXGBoost Results:\n")
print(xgb_cm)
cat("XGBoost AUC:", round(xgb_auc, 4), "\n")

# Save XGBoost importance
xgb_importance <- xgb.importance(
  feature_names = colnames(X_train),
  model = xgb_model
)

xgb_importance_annotated <- xgb_importance %>%
  left_join(
    selected_feature_table,
    by = c("Feature" = "ID")
  )

write.csv(
  xgb_importance_annotated,
  "results/XGBoost_feature_importance.csv",
  row.names = FALSE
)

# ============================================================
# 9. Model comparison summary
# ============================================================

model_summary <- data.frame(
  Model = c("LASSO", "Random Forest", "XGBoost"),
  Accuracy = c(
    lasso_cm$overall["Accuracy"],
    rf_cm$overall["Accuracy"],
    xgb_cm$overall["Accuracy"]
  ),
  Sensitivity = c(
    lasso_cm$byClass["Sensitivity"],
    rf_cm$byClass["Sensitivity"],
    xgb_cm$byClass["Sensitivity"]
  ),
  Specificity = c(
    lasso_cm$byClass["Specificity"],
    rf_cm$byClass["Specificity"],
    xgb_cm$byClass["Specificity"]
  ),
  AUC = c(lasso_auc, rf_auc, xgb_auc)
)

write.csv(
  model_summary,
  "results/ML_model_comparison_summary.csv",
  row.names = FALSE
)

cat("\nModel comparison summary:\n")
print(model_summary)

# ============================================================
# 10. ROC curve plot
# ============================================================

roc_df <- rbind(
  data.frame(
    FPR = 1 - lasso_roc$specificities,
    TPR = lasso_roc$sensitivities,
    Model = paste0("LASSO (AUC = ", round(lasso_auc, 3), ")")
  ),
  data.frame(
    FPR = 1 - rf_roc$specificities,
    TPR = rf_roc$sensitivities,
    Model = paste0("Random Forest (AUC = ", round(rf_auc, 3), ")")
  ),
  data.frame(
    FPR = 1 - xgb_roc$specificities,
    TPR = xgb_roc$sensitivities,
    Model = paste0("XGBoost (AUC = ", round(xgb_auc, 3), ")")
  )
)

roc_plot <- ggplot(roc_df, aes(x = FPR, y = TPR, color = Model)) +
  geom_line(linewidth = 1) +
  geom_abline(linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "ROC Curves for MDD Classification Models",
    x = "False Positive Rate",
    y = "True Positive Rate"
  )

ggsave(
  "plots/ML_ROC_curve.png",
  roc_plot,
  width = 7,
  height = 5,
  dpi = 300
)

# ============================================================
# 11. Top feature importance plot
# ============================================================

top_rf_features <- rf_importance %>%
  filter(!is.na(`Gene Symbol`), `Gene Symbol` != "") %>%
  head(20)

rf_importance_plot <- ggplot(
  top_rf_features,
  aes(x = reorder(`Gene Symbol`, Importance), y = Importance)
) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Random Forest Top 20 Feature Importance",
    x = "Gene Symbol",
    y = "Importance"
  )

ggsave(
  "plots/RandomForest_top20_feature_importance.png",
  rf_importance_plot,
  width = 7,
  height = 5,
  dpi = 300
)

cat("\nMachine learning model comparison completed successfully.\n")