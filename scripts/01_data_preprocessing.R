# ============================================================
# MDD Gene Expression Analysis — GSE98793
# 01: Data Loading & Preprocessing
# ============================================================

library(GEOquery)
library(Biobase)
library(dplyr)

# ---- 1. Load GEO series matrix ----
gset <- getGEO(filename = "data/raw/GSE98793_series_matrix.txt")

# If multiple platforms are returned, use the first ExpressionSet
if (is.list(gset)) {
  gset <- gset[[1]]
}

# ---- 2. Extract expression matrix, phenotype data, and feature data ----
expr_mat <- exprs(gset)
pheno_data <- pData(gset)
feature_data <- fData(gset)

cat("Expression matrix dimensions (features x samples):", dim(expr_mat), "\n")
cat("\nPhenotype columns:\n")
print(colnames(pheno_data))

cat("\nFeature annotation columns:\n")
print(colnames(feature_data))

# ---- 3. Inspect phenotype fields for group labels ----
cat("\ncharacteristics_ch1 values:\n")
print(table(pheno_data$characteristics_ch1, useNA = "ifany"))

if ("characteristics_ch1.1" %in% colnames(pheno_data)) {
  cat("\ncharacteristics_ch1.1 values:\n")
  print(table(pheno_data$characteristics_ch1.1, useNA = "ifany"))
}

# ---- 4. Build group labels safely ----
group <- case_when(
  grepl("case", pheno_data$characteristics_ch1, ignore.case = TRUE) ~ "MDD",
  grepl("control|healthy", pheno_data$characteristics_ch1, ignore.case = TRUE) ~ "HC",
  TRUE ~ NA_character_
)

if (any(is.na(group)) && "characteristics_ch1.1" %in% colnames(pheno_data)) {
  idx <- which(is.na(group))
  group[idx] <- case_when(
    grepl("case", pheno_data$characteristics_ch1.1[idx], ignore.case = TRUE) ~ "MDD",
    grepl("control|healthy", pheno_data$characteristics_ch1.1[idx], ignore.case = TRUE) ~ "HC",
    TRUE ~ NA_character_
  )
}

pheno_factor <- factor(group, levels = c("HC", "MDD"))

cat("\nGroup counts:\n")
print(table(pheno_factor, useNA = "ifany"))

if (any(is.na(pheno_factor))) {
  stop("Some samples could not be assigned to HC or MDD. Check phenotype fields before continuing.")
}

# ---- 5. Check expression range and log2 transform only if needed ----
expr_range <- range(expr_mat, na.rm = TRUE)
cat("\nExpression range:", expr_range, "\n")

if (max(expr_mat, na.rm = TRUE) > 100) {
  expr_mat <- log2(expr_mat + 1)
  cat("Applied log2 transformation\n")
} else {
  cat("Data appears to already be on a log2 scale; no additional log2 transformation applied\n")
}

# ---- 6. Remove rows with missing values ----
cat("\nTotal missing values:", sum(is.na(expr_mat)), "\n")
expr_mat <- expr_mat[complete.cases(expr_mat), ]
cat("Dimensions after removing rows with NA:", dim(expr_mat), "\n")

# ---- 7. Filter low-signal features ----
probe_mean <- rowMeans(expr_mat)
keep <- probe_mean > 5
expr_filtered <- expr_mat[keep, ]

cat("\nFeatures before filtering:", nrow(expr_mat), "\n")
cat("Features after filtering:", nrow(expr_filtered), "\n")

# ---- 8. Sanity checks for sample alignment ----
stopifnot(length(pheno_factor) == ncol(expr_filtered))

if ("geo_accession" %in% colnames(pheno_data)) {
  stopifnot(all(colnames(expr_filtered) == pheno_data$geo_accession))
}

# ---- 9. Create design matrix for limma ----
design <- model.matrix(~ pheno_factor)

cat("\nDesign matrix preview:\n")
print(head(design))

# ---- 10. Save preprocessed objects ----
save(
  expr_filtered,
  pheno_factor,
  pheno_data,
  feature_data,
  design,
  file = "data/processed/MDD_preprocessed.RData"
)

cat("\nSaved: data/processed/MDD_preprocessed.RData\n")
