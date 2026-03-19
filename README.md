# MDD Gene Expression & GWAS Analysis
Machine learning analysis of Major Depressive Disorder using RNA-seq and GWAS data(R)

Computational investigation of Major Depressive Disorder (MDD) using 
RNA-seq gene expression data and genome-wide SNP data. The analysis 
pipeline covers end-to-end preprocessing, statistical modeling, network 
analysis, and machine learning classification.

---

## Background

Major Depressive Disorder is a complex psychiatric condition with significant 
genetic components. This project applies bioinformatics and machine learning 
methods to real patient-level data to identify molecular signatures associated 
with MDD — including differentially expressed genes, co-expressed gene modules, 
and disease-associated SNPs.

---

## Pipeline

### 1. Data Preprocessing
- RNA-seq expression data: 8,923 genes across MDD and healthy control (HC) subjects
- Applied quantile normalization and log2 transformation for cross-sample comparability
- Filtered low-signal genes using coefficient of variation (CoV) threshold

### 2. Differential Expression & Feature Selection
- Logistic regression across all genes to rank MDD-associated candidates
- LASSO regression (glmnet) for sparse, high-dimensional feature selection
- **Classification accuracy: 81.53%**

### 3. ML Model Benchmarking
Compared three classifiers on the MDD vs. HC task:

| Model | Performance |
|-------|-------------|
| LASSO (baseline) | 81.53% |
| Random Forest | ↑ higher |
| XGBoost | ↑ higher |

Random Forest and XGBoost both outperformed the LASSO baseline.

### 4. Gene Co-expression Network
- Constructed gene network using Pearson correlation with Random Matrix Theory threshold
- Community detection via Louvain/Leiden algorithm to identify functional gene modules
- Network visualization with igraph

### 5. GWAS
- Tested SNP–phenotype associations using Fisher's exact test and logistic regression
- Manhattan plot generated to visualize genome-wide significance patterns

---

## Results Summary

- Identified a subset of genes significantly differentially expressed in MDD vs. HC
- LASSO selected a sparse gene signature with 81.53% classification accuracy
- Random Forest and XGBoost further improved classification performance
- Co-expression network revealed two distinct gene modules (clust1, clust2)
- GWAS flagged candidate SNPs with significant phenotype association

---

## Libraries

`glmnet` · `ranger` · `xgboost` · `preprocessCore` · `igraph` · 
`snpStats` · `ggplot2` · `dplyr` · `broom`

---

## File Structure
```
├── 01_data_preprocessing.R
├── 02_differential_expression.R
├── 03_feature_selection_LASSO.R
├── 04_coexpression_network.R
├── 05_GWAS_analysis.R
├── 06_model_comparison.R
└── README.md
```
