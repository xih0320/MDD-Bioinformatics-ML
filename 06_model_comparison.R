# ============================================================
# MDD Gene Expression Analysis
# 06: ML Model Comparison — LASSO vs Random Forest vs XGBoost
# ============================================================
# Requires: GxS_covfilter, pheno_factor from 01_data_preprocessing.R

library(glmnet)
library(ranger)
library(xgboost)
library(caret)

# Prepare matrix input (subjects as rows)
X = t(GxS_covfilter)
y = as.numeric(pheno_factor) - 1  # 0 = HC, 1 = MDD

set.seed(42)

# ---- 1. LASSO (baseline) ----
glmnet_model = cv.glmnet(x = X, y = y, alpha = 1,
                         family = "binomial", type.measure = "class", keep = TRUE)

lasso_predicted = predict(glmnet_model, s = glmnet_model$lambda.min,
                           type = "class", newx = X)
lasso_accuracy  = mean(as.integer(lasso_predicted) == y)

cat("LASSO Classification Accuracy:", round(lasso_accuracy * 100, 2), "%\n")
cat("Confusion matrix:\n")
print(table(Predicted = lasso_predicted, Actual = y))

# ---- 2. Random Forest ----
rf_fit = ranger(
  y ~ .,
  data           = data.frame(X, y = as.factor(y)),
  num.trees      = 1000,
  classification = TRUE,
  importance     = "permutation",
  seed           = 42
)

rf_accuracy = 1 - rf_fit$prediction.error
cat("\nRandom Forest OOB Accuracy:", round(rf_accuracy * 100, 2), "%\n")
cat("Confusion matrix:\n")
print(rf_fit$confusion.matrix)

# Variable importance (top 20 genes)
rf_imp = sort(rf_fit$variable.importance, decreasing = TRUE)
barplot(rf_imp[1:20], las = 2, main = "Random Forest — Top 20 Gene Importance",
        col = "steelblue", cex.names = 0.7)

# ---- 3. XGBoost ----
dtrain = xgb.DMatrix(data = X, label = y)

params = list(
  objective  = "binary:logistic",
  eval_metric = "error",
  booster    = "gbtree",
  eta        = 0.1,
  max_depth  = 4
)

xgb_model = xgboost(params = params, data = dtrain, nrounds = 50,
                     verbose = 0)

xgb_probs      = predict(xgb_model, dtrain)
xgb_predicted  = as.integer(xgb_probs >= 0.5)
xgb_accuracy   = mean(xgb_predicted == y)

cat("\nXGBoost Training Accuracy:", round(xgb_accuracy * 100, 2), "%\n")
cat("Confusion matrix:\n")
print(table(Predicted = xgb_predicted, Actual = y))

# ---- Summary ----
cat("\n========== Model Comparison ==========\n")
cat(sprintf("%-20s %.2f%%\n", "LASSO",         lasso_accuracy * 100))
cat(sprintf("%-20s %.2f%%\n", "Random Forest", rf_accuracy * 100))
cat(sprintf("%-20s %.2f%%\n", "XGBoost",       xgb_accuracy * 100))
