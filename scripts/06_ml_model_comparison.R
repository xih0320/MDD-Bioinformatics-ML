# ============================================================
# 06: Machine Learning Model Comparison (FINAL FIXED VERSION)
# LASSO vs Random Forest vs XGBoost + Gene Symbol FIX
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

# ---- 2. Load data ----
expr_matrix <- readRDS("results/batch_corrected_expression_matrix.rds")
pheno_factor <- readRDS("results/pheno_factor.rds")

de_results <- read.csv("results/MDD_DE_results_full.csv", check.names = FALSE)

# ---- 3. Feature selection ----
sig_features <- de_results %>%
  filter(adj.P.Val < 0.05) %>%
  pull(ID) %>%
  unique()

cat("Significant DE genes:", length(sig_features), "\n")

common_features <- intersect(sig_features, rownames(expr_matrix))

# fallback
if(length(common_features) < 10){
  top_features <- de_results %>%
    arrange(adj.P.Val) %>%
    pull(ID) %>%
    unique() %>%
    head(100)
  
  common_features <- intersect(top_features, rownames(expr_matrix))
}

cat("Final feature count:", length(common_features), "\n")


selected_feature_table <- de_results %>%
  filter(ID %in% common_features) %>%
  select(ID, `Gene Symbol`, logFC, adj.P.Val)

# ---- 4. ML matrix ----
X <- t(expr_matrix[common_features, ])
y <- factor(pheno_factor, levels = c("HC", "MDD"))

# ---- 5. Split ----
set.seed(42)
train_index <- createDataPartition(y, p = 0.8, list = FALSE)

X_train <- X[train_index, ]
X_test  <- X[-train_index, ]

y_train <- y[train_index]
y_test  <- y[-train_index]

y_train_num <- ifelse(y_train == "MDD", 1, 0)
y_test_num  <- ifelse(y_test == "MDD", 1, 0)

# ============================================================
# 6. LASSO
# ============================================================

lasso_cv <- cv.glmnet(
  x = X_train,
  y = y_train_num,
  family = "binomial",
  alpha = 1,
  type.measure = "auc"
)

lasso_prob <- predict(lasso_cv, newx = X_test, s = "lambda.min", type = "response") %>% as.numeric()

lasso_roc <- roc(y_test_num, lasso_prob, quiet = TRUE)
lasso_auc <- auc(lasso_roc)

# ============================================================
# 7. Random Forest
# ============================================================

rf_train_df <- as.data.frame(X_train, check.names = FALSE)
rf_train_df$Group <- y_train

rf_test_df <- as.data.frame(X_test, check.names = FALSE)

set.seed(42)

rf_model <- ranger(
  dependent.variable.name = "Group",
  data = rf_train_df,
  num.trees = 500,
  probability = TRUE,
  importance = "permutation",
  seed = 42
)

rf_prob <- predict(rf_model, data = rf_test_df)$predictions[, "MDD"]
rf_roc <- roc(y_test_num, rf_prob, quiet = TRUE)
rf_auc <- auc(rf_roc)

# ---- FIX: Gene Symbol join ----
rf_importance <- data.frame(
  Feature = names(rf_model$variable.importance),
  Importance = as.numeric(rf_model$variable.importance),
  check.names = FALSE
)

rf_importance$Feature <- as.character(rf_importance$Feature)
selected_feature_table$ID <- as.character(selected_feature_table$ID)

rf_importance <- rf_importance %>%
  left_join(selected_feature_table, by = c("Feature" = "ID")) %>%
  mutate(GeneLabel = ifelse(
    !is.na(`Gene Symbol`) & `Gene Symbol` != "",
    `Gene Symbol`,
    Feature
  )) %>%
  arrange(desc(Importance))

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

xgb_model <- xgb.train(
  params = list(objective = "binary:logistic", eval_metric = "auc"),
  data = dtrain,
  nrounds = 100,
  verbose = 0
)

xgb_prob <- predict(xgb_model, dtest)
xgb_roc <- roc(y_test_num, xgb_prob, quite = TRUE)
xgb_auc <- auc(xgb_roc)

# ============================================================
# 9. ROC plot
# ============================================================

roc_df <- rbind(
  data.frame(FPR = 1 - lasso_roc$specificities, TPR = lasso_roc$sensitivities, Model = "LASSO"),
  data.frame(FPR = 1 - rf_roc$specificities, TPR = rf_roc$sensitivities, Model = "Random Forest"),
  data.frame(FPR = 1 - xgb_roc$specificities, TPR = xgb_roc$sensitivities, Model = "XGBoost")
)

roc_plot <- ggplot(roc_df, aes(FPR, TPR, color = Model)) +
  geom_line(size = 1.2) +
  geom_abline(linetype = "dashed") +
  theme_minimal() +
  labs(
    title = paste0(
      "ROC Curves (LASSO=", round(lasso_auc,2),
      ", RF=", round(rf_auc,2),
      ", XGB=", round(xgb_auc,2), ")"
    )
  )

ggsave("plots/ML_ROC_curve.png", roc_plot, width = 7, height = 5)

# ============================================================
# 10. Feature importance plot
# ============================================================

top_rf <- rf_importance %>% head(20)

rf_plot <- ggplot(top_rf, aes(x = reorder(GeneLabel, Importance), y = Importance)) +
  geom_col(fill = "#2C7BB6") +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Top 20 Important Genes (Random Forest)",
    x = "Gene Symbol",
    y = "Importance"
  )

ggsave("plots/RandomForest_top20_feature_importance.png", rf_plot, width = 7, height = 5)

