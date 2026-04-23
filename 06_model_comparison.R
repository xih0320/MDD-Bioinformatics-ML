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

# ---- Train/Test Split (80/20) ----
set.seed(42)
train_idx = createDataPartition(y, p = 0.8, list = FALSE)
X_train = X[train_idx, ]
X_test  = X[-train_idx, ]
y_train = y[train_idx]
y_test  = y[-train_idx]

cat("Train size:", length(y_train), "| Test size:", length(y_test), "\n")

# ---- 1. LASSO (baseline) ----
set.seed(42)
glmnet_model = cv.glmnet(x = X_train, y = y_train, alpha = 1,
                         family = "binomial", type.measure = "class", keep = TRUE)

lasso_predicted = predict(glmnet_model, s = glmnet_model$lambda.min,
                          type = "class", newx = X_test)
lasso_accuracy  = mean(as.integer(lasso_predicted) == y_test)

cat("\nLASSO Test Accuracy:", round(lasso_accuracy * 100, 2), "%\n")
cat("Confusion matrix:\n")
print(table(Predicted = lasso_predicted, Actual = y_test))

# ---- 2. Random Forest ----
set.seed(42)
rf_fit = ranger(
  y ~ .,
  data           = data.frame(X_train, y = as.factor(y_train)),
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
dtrain = xgb.DMatrix(data = X_train, label = y_train)
dtest  = xgb.DMatrix(data = X_test,  label = y_test)

set.seed(42)
xgb_model = xgb.train(
  params = list(
    objective    = "binary:logistic",
    eval_metric  = "error",
    learning_rate = 0.1,
    max_depth    = 4
  ),
  data    = dtrain,
  nrounds = 50
)

xgb_probs     = predict(xgb_model, dtest)
xgb_predicted = as.integer(xgb_probs >= 0.5)
xgb_accuracy  = mean(xgb_predicted == y_test)

cat("\nXGBoost Test Accuracy:", round(xgb_accuracy * 100, 2), "%\n")
cat("Confusion matrix:\n")
print(table(Predicted = xgb_predicted, Actual = y_test))

# ---- Summary ----
cat("\n========== Model Comparison (Test Set) ==========\n")
cat(sprintf("%-20s %.2f%%\n", "LASSO",         lasso_accuracy * 100))
cat(sprintf("%-20s %.2f%%\n", "Random Forest (OOB)", rf_accuracy * 100))
cat(sprintf("%-20s %.2f%%\n", "XGBoost",       xgb_accuracy * 100))
