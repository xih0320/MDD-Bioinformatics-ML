# MDD Gene Expression & GWAS Analysis Report

## Dataset

- RNA-seq expression data: 8,923 genes, MDD patients vs healthy controls (HC)
- GWAS data: 89 subjects, 34 SNPs (PLINK format)

---

## 1. Data Preprocessing

- Quantile normalization + log2 transformation applied across all samples
- Coefficient of variation (CoV) filter removed low-signal genes
- Result: 8,923 genes filtered down to 5,587 genes

---

## 2. Differential Expression

Welch t-test applied across all 5,587 genes. Top results:

| Gene | t-statistic | p-value |
|------|-------------|---------|
| MDGA1 | -4.603 | 8.77e-06 |
| ZDHHC20 | -4.147 | 5.56e-05 |
| NPFF | 3.921 | 1.32e-04 |
| ARFGAP1 | -3.849 | 1.75e-04 |
| FAM138A | 3.816 | 1.98e-04 |

Top gene **MDGA1** is involved in synaptic organization and relevant to psychiatric disorders.

---

## 3. Pathway Enrichment (MSigDB Reactome)

Top 200 DE genes submitted to MSigDB Reactome:

| Pathway | Overlap (k) | p-value | FDR |
|---------|-------------|---------|-----|
| Post-translational protein modification | 23 | 4.53e-7 | 7.86e-4 |
| Metabolism of RNA | 15 | 1.77e-6 | 1.54e-3 |
| SARS-CoV-2 Infection | 9 | 1.24e-5 | 7.18e-3 |
| Infectious disease | 16 | 5.15e-5 | 1.79e-2 |

---

## 4. Gene Co-expression Network

- Pearson correlation matrix built across filtered genes
- Random Matrix Theory threshold applied to remove noise edges
- Community detection (Fast-greedy, Louvain, Leiden) identified gene modules
- Two major modules exported for pathway enrichment
- UMAP visualization showed continuous biological structure across gene modules

---

## 5. GWAS Results

Fisher exact test and logistic regression across 34 SNPs:

- Top SNP: **rs7835221** (chromosome 8)
- Strong genotype-phenotype association confirmed by logistic regression
- rs7835221 exceeded genome-wide significance threshold in Manhattan plot

---

## 6. ML Model Comparison

| Model | Accuracy |
|-------|----------|
| LASSO (cv.glmnet) | 81.53% |
| Random Forest | higher than LASSO |
| XGBoost | higher than LASSO |

LASSO confusion matrix:

|  | Predicted HC | Predicted MDD |
|--|-------------|---------------|
| Actual HC | 62 | 17 |
| Actual MDD | 12 | 66 |

---

## Summary

| Analysis | Result |
|----------|--------|
| Genes retained after filtering | 5,587 / 8,923 |
| Top DE gene | MDGA1 (p = 8.77e-06) |
| Top enriched pathway | Post-translational protein modification |
| LASSO accuracy | 81.53% |
| Top GWAS SNP | rs7835221 (chromosome 8) |
