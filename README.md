# Gene Expression Analysis and Machine Learning Modeling for Major Depressive Disorder (MDD)

An end-to-end bioinformatics project investigating gene expression changes associated with Major Depressive Disorder (MDD) using public GEO microarray data (GSE98793, n = 192).

---

## Overview

This project integrates statistical analysis, functional enrichment, protein-protein interaction (PPI) network analysis, and machine learning to characterize transcriptomic changes in MDD.

The workflow spans from raw expression data to biologically interpretable findings, with a focus on identifying immune-related gene signatures and evaluating their predictive power.

---

## Dataset

* **GEO Accession:** GSE98793
* **Samples:** 192 (64 healthy controls, 128 MDD cases)
* **Platform:** Affymetrix microarray
* **Source:** Peripheral blood gene expression data

---

## Pipeline

### 1. Data Preprocessing (01_data_preprocessing.R)

* Loaded GEO series matrix using GEOquery
* Applied log2 transformation
* Filtered low-expression probes (mean > 5)
* Extracted phenotype and annotation data

### 2. QC and Batch Correction (02_qc_pca.R)

* PCA before/after batch correction
* Batch effect removed using ComBat (sva)
* Reduced PC1 variance from 29.7% → 20.0%

### 3. Differential Expression (03_differential_expression.R)

* limma with empirical Bayes
* Model: `~ pheno_factor + batch`
* Threshold: adj.P.Val < 0.05, |logFC| > 0.3
* Identified **26 significant DEGs**

### 4. Functional Enrichment (04_GO_KEGG_enrichment.R)

* GO enrichment (clusterProfiler)
* KEGG pathway analysis
* Converted gene symbols → Entrez IDs

### 5. PPI Network (05_string_ppi_network.R)

* STRING v12 (score ≥ 400)
* Network constructed via igraph
* Hub genes: ARG1, LCN2, HP, S100A12

### 6. Machine Learning (06_ml_model_comparison.R)

* Selected **150 significant genes** as features
* Models: LASSO, Random Forest, XGBoost
* Train/test split: 80/20
* Evaluation: Accuracy, AUC, Sensitivity, Specificity

---

## Key Results

### Differential Expression

* 26 DEGs identified
* Top genes: RORA, GZMK, RETN, MAFG

### Enrichment Analysis

* Significant GO terms enriched in:

  * vesicle lumen
  * secretory granule lumen
* Indicates **neutrophil-mediated immune activity**

### PPI Network

* Hub genes:

  * ARG1 (degree = 8)
  * LCN2 (degree = 8)
* Suggests coordinated innate immune activation

---

## Machine Learning Performance

* **LASSO AUC: 0.92 (best)**
* Random Forest AUC: 0.84
* XGBoost AUC: 0.80

LASSO outperformed other models, likely due to its ability to handle high-dimensional data and perform embedded feature selection.

---

## Biological Interpretation

Top features identified by Random Forest include:

* CD48
* SH2D1A
* GZMK
* RETN

These genes are involved in **immune cell activation, inflammation, and T-cell signaling**.

Importantly, these findings are consistent with:

* GO enrichment results (immune-related pathways)
* PPI network (innate immune hub genes)

> Together, both statistical and machine learning analyses independently highlight **immune dysregulation as a key mechanism in MDD**.

---

## Project Structure

```
MDD_microarray_analysis/
├── scripts/
├── results/
├── plots/
```

---

## Reproducibility

```r
setwd("your_project_directory")

source("scripts/01_data_preprocessing.R")
source("scripts/02_qc_pca.R")
source("scripts/03_differential_expression.R")
source("scripts/04_GO_KEGG_enrichment.R")
source("scripts/05_string_ppi_network.R")
source("scripts/06_ml_model_comparison.R")
```

---

## Tools & Packages

* GEOquery, limma, sva
* clusterProfiler, org.Hs.eg.db
* STRINGdb, igraph
* ranger, glmnet, xgboost
* ggplot2

---

## Conclusion

This project demonstrates an integrated bioinformatics and machine learning workflow to identify biologically meaningful patterns in MDD.

Both statistical and predictive modeling approaches consistently indicate that **immune-related mechanisms play a central role in MDD**, highlighting potential biomarkers for future research.

---

## Author

**Xinyue (Kitty) Hu**
M.S. Computer Science, University of Tulsa
GitHub: github.com/xih0320
