# Gene Expression Analysis and Machine Learning Modeling for Major Depressive Disorder (MDD)

An end-to-end transcriptomic analysis project investigating gene expression changes associated with Major Depressive Disorder(MDD) using public GEO data(GSE98793, n=192).

---

## Overview

This project integrates statistical modeling, functional enrichment, and protein-protein interaction (PPI) network analysis with machine learning approaches (LASSO, Random Forest, XGBoost) to predict disease status from gene expression profiles to characterize transcriptomic changes in MDD. The analysis spans from raw expression data to biologically interpretable findings, with a focus on immune-related gene signatures.

---

## Dataset

- **GEO Accession:** GSE98793
- **Samples:** 192 (64 healthy controls, 128 MDD cases)
- **Platform:** Affymetrix microarray
- **Source:** Blood-based gene expression data

---

## Pipeline

### 1. Data Preprocessing (`01_data_preprocessing.R`)
- Loaded series matrix from GEO using `GEOquery`
- Applied log2 transformation where needed
- Filtered low-signal probes (mean expression > 5)
- Extracted phenotype and feature annotation data

### 2. QC and Batch Correction (`02_qc_pca.R`)
- PCA before and after batch correction, colored by phenotype and batch
- Batch effect removed using `ComBat` (sva package)
- Batch correction reduced PC1 variance from 29.7% to 20.0%, confirming successful removal

### 3. Differential Expression Analysis (`03_differential_expression.R`)
- Used `limma` with empirical Bayes moderation
- Design matrix included batch as a covariate: `~ pheno_factor + batch`
- BH-adjusted p-value < 0.05, |logFC| > 0.3
- **26 significant DEGs identified**
- Volcano plot generated for visualization

### 4. Functional Enrichment (`04_GO_KEGG_enrichment.R`)
- Gene Ontology enrichment using `clusterProfiler` (`ont = "ALL"`)
- KEGG pathway analysis
- Gene symbols converted to Entrez IDs via `bitr`

### 5. PPI Network Analysis (`05_string_ppi_network.R`)
- Queried STRING database (v12.0, score threshold = 400) via `STRINGdb`
- Network constructed with `igraph`
- Hub genes identified by degree centrality

### 6. Machine Learning Modeling (`06_ml_model_comparison.R`)
- Selected top 200 DE genes as features
- Built classification models: LASSO (glmnet), Random Forest (ranger), XGBoost
- Data split into training (80%) and test (20%) sets
- Model performance evaluated using accuracy, sensitivity, specificity, and AUC
- LASSO achieved the best performance (AUC = 0.94), demonstrating strong predictive power and feature sparsity

---

## Results

**Differential Expression**
- 26 significant DEGs under strict thresholds (adj.P.Val < 0.05, |logFC| > 0.3)
- Top genes include RORA, GZMK, RETN, MAFG
- A relatively relaxed logFC threshold(|logFC| > 0.3) was used to retain sufficient genes for downstream enrichment and network analysis, given the moderate effect sizes typically observed in blood-based transcriptomic data.

**GO Enrichment**
- 6 significant GO terms identified, all under Cellular Component (CC)
- Enriched terms: vesicle lumen, secretory granule lumen, cytoplasmic vesicle lumen
- Core genes: HP, LCN2, ARG1, OLFM4
- Suggests dysregulation of neutrophil-mediated immune secretion in MDD

**PPI Network**
- Hub genes by degree: ARG1 (8), LCN2 (8), HP (6), S100A12 (6)
- Network topology consistent with coordinated innate immune activation

**Biological Interpretation**

The most consistently enriched signal points to immune dysregulation — specifically neutrophil granule components and vesicle-mediated secretion. This aligns with existing literature implicating peripheral inflammatory processes in MDD pathophysiology. RORA, a circadian rhythm regulator with known associations to depression, also appeared among top DEGs. While the number of DEGs is relatively small, the consistency across enrichment and network analyses strengthens the confidence in these immune-related signals.

**Machine Learning Performance**
- LASSO achieved the best performance (Accuracy = 0.865, AUC = 0.94)
- Random Forest and XGBoost showed lower specificity and overall AUC
- Results suggest that sparse linear models are particularly effective for high-dimensional gene expression data

---

## Project Structure

```
MDD_RNAseq_analysis/
├── scripts/
│   ├── 01_data_preprocessing.R
│   ├── 02_qc_pca.R
│   ├── 03_differential_expression.R
│   ├── 04_GO_KEGG_enrichment.R
│   ├── 05_string_ppi_network.R
│   └── 06_ml_model_comparison.R

├── results/
│   ├── MDD_DE_genes_sig.csv
│   ├── MDD_DE_results_full.csv
│   ├── GO_enrichment_results.csv
│   ├── KEGG_enrichment_results.csv
│   ├── STRING_hub_genes.csv
│   │
│   ├── ML_model_comparison_summary.csv
│   ├── LASSO_selected_features.csv
│   ├── XGBoost_feature_importance.csv
│   └── RandomForest_feature_importance.csv

├── plots/
│   ├── MDD_volcano_plot.png
│   ├── GO_dotplot.png
│   ├── STRING_PPI_network.png
│   ├── STRING_hub_genes_barplot.png
│   │
│   ├── ML_ROC_curve.png
│   └── RandomForest_top20_feature_importance.png
```

---

## Reproducibility

All scripts use relative paths. To reproduce the full pipeline:

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

| Category | Tools |
|---|---|
| Data retrieval | GEOquery, Biobase |
| Batch correction | sva (ComBat) |
| Differential expression | limma |
| Enrichment analysis | clusterProfiler, org.Hs.eg.db, enrichplot |
| Network analysis | STRINGdb, igraph |
| Visualization | ggplot2 |

---

## Future Directions

- Cross-dataset validation using independent MDD cohorts from GEO
- GSEA to capture pathway-level signals beyond individual DEGs
- Incorporation of clinical covariates (e.g., medication status, symptom severity) into the model

---

## Author

Xinyue (Kitty) Hu  
M.S. Computer Science, University of Tulsa  
GitHub: [github.com/xih0320](https://github.com/xih0320)
