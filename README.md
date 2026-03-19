# MDD Gene Expression & GWAS Analysis

## Highlights
- Built end-to-end bioinformatics pipeline integrating RNA-seq and GWAS data  
- Applied LASSO, Random Forest, and XGBoost for high-dimensional classification  
- Identified gene signatures and SNP associations linked to Major Depressive Disorder  

---

## Overview
This project performs a computational analysis of Major Depressive Disorder (MDD) using RNA-seq gene expression data and genome-wide SNP data.

The pipeline integrates statistical modeling, network analysis, and machine learning to identify molecular signatures associated with MDD and evaluate predictive performance.

---

## Why it matters
Major Depressive Disorder is a complex psychiatric condition with strong genetic components.  

This project demonstrates how machine learning and bioinformatics methods can be used to:
- Identify disease-associated gene expression patterns  
- Discover SNP-level associations  
- Build predictive models for disease classification  

---

## Pipeline

### 1. Data Preprocessing
- RNA-seq dataset: 8,923 genes across MDD and healthy control (HC) subjects  
- Quantile normalization + log2 transformation  
- Filtered low-variance genes using coefficient of variation (CoV)  

---

### 2. Differential Expression & Feature Selection
- Logistic regression applied across all genes to rank MDD-associated candidates  
- LASSO regression (glmnet) for sparse feature selection in high-dimensional space  
- Achieved **81.53% classification accuracy**  

---

### 3. Model Benchmarking
Compared multiple machine learning models for MDD classification:

| Model | Performance |
|-------|------------|
| LASSO | 81.53% |
| Random Forest | Higher |
| XGBoost | Higher |

Tree-based models outperformed linear baselines, indicating non-linear gene interactions.

---

### 4. Gene Co-expression Network
- Constructed gene network using Pearson correlation + Random Matrix Theory threshold  
- Applied Louvain/Leiden algorithm for community detection  
- Identified two major functional gene modules  
- Visualization using igraph  

---

### 5. GWAS Analysis
- SNP–phenotype association testing using Fisher’s exact test and logistic regression  
- Generated Manhattan plot for genome-wide significance visualization  

---

## Results Summary
- Identified differentially expressed genes associated with MDD  
- LASSO selected a compact gene signature with strong predictive performance  
- Random Forest and XGBoost improved classification accuracy  
- Co-expression analysis revealed distinct gene modules  
- GWAS identified candidate SNPs associated with MDD  

---

## Tech Stack
R · glmnet · ranger · xgboost · preprocessCore · igraph · snpStats · ggplot2 · dplyr  

---

## File Structure
├── 01_data_preprocessing.R
├── 02_differential_expression.R
├── 03_feature_selection_LASSO.R
├── 04_coexpression_network.R
├── 05_GWAS_analysis.R
├── 06_model_comparison.R
└── README.md
