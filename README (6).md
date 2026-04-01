# MDD Gene Expression & GWAS Analysis

## Highlights
- Built end-to-end bioinformatics pipeline integrating RNA-seq and GWAS data
- Applied LASSO, Random Forest, and XGBoost for high-dimensional classification
- Identified gene signatures and SNP associations linked to Major Depressive Disorder

---

## Overview
This project performs a computational analysis of Major Depressive Disorder (MDD) using RNA-seq gene expression data and genome-wide SNP data from 157 matched subjects (80 MDD patients and 80 healthy controls).

The pipeline integrates statistical modeling, network analysis, and machine learning to identify molecular signatures associated with MDD and evaluate predictive performance.

---

## Why it matters
Major Depressive Disorder is a complex psychiatric condition with strong genetic components.

This project demonstrates how machine learning and bioinformatics methods can be used to:
- Identify disease-associated gene expression patterns
- Discover SNP-level associations
- Build predictive models for disease classification

---

## Data
RNA-seq expression data and clinical phenotype data from GEO dataset GSE98793.
Data is not included in this repository. Download from:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98793

---

## Pipeline

### 1. Data Preprocessing
- RNA-seq dataset: 8,923 genes across 157 matched subjects (80 MDD, 80 HC)
- Quantile normalization to remove batch effects between samples
- Log2 transformation to stabilize variance
- Filtered genes using coefficient of variation (CoV < 0.045), retaining 5,587 genes

---

### 2. Differential Expression & Feature Selection
- t-test applied across all 5,587 genes to identify MDD-associated candidates
- Top differentially expressed gene: **MDGA1** (p = 8.77e-06)
- No genes passed FDR < 0.05 threshold, motivating the use of top variable genes as ML features
- LASSO regression (glmnet) for sparse feature selection, selecting 16 key genes
- LASSO in-sample classification accuracy: **~80%**

---

### 3. Model Benchmarking
Compared multiple machine learning models for MDD classification using 80/20 train/test split:

| Model | Accuracy | AUC | Evaluation Method |
|--------|----------|-----|-------------------|
| LASSO | — | 0.679 | In-sample |
| Random Forest | 59.52% | 0.679 | OOB (Out-of-Bag) |
| XGBoost | **67.74%** | **0.726** Best | Test set |

Note: Low accuracy reflects the challenge of high-dimensional classification with small sample size (n=157).

Note: Low accuracy reflects the challenge of high-dimensional classification with small sample size (n=157). LASSO in-sample accuracy is optimistic due to no held-out test set.

---

### 4. Gene Co-expression Network
- Constructed gene co-expression network using Pearson correlation
- Applied Random Matrix Theory (RMT) to determine correlation threshold (0.45)
- Applied Fast-Greedy algorithm for community detection
- Identified two major functional gene modules (clust1.txt, clust2.txt)
- Visualization using igraph

---

### 5. GWAS Analysis
- SNP–phenotype association testing using logistic regression
- Top associated SNP: **rs7835221** (p = 0.025)
- Generated Manhattan plot for genome-wide significance visualization

---

## Results Summary
- Top differentially expressed gene: MDGA1 (p = 8.77e-06)
- LASSO selected 16 key gene biomarkers including MDGA1, OSBPL3, PANX1, NPFF
- Co-expression analysis revealed two distinct gene modules for pathway enrichment
- GWAS identified rs7835221 as top candidate SNP associated with MDD

---

## Tech Stack
R · glmnet · ranger · xgboost · preprocessCore · igraph · snpStats · ggplot2 · dplyr

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
