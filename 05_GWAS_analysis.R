# ============================================================
# MDD Gene Expression Analysis
# 05: GWAS — SNP-Phenotype Association Analysis
# ============================================================

library(snpStats)
library(dplyr)
library(ggplot2)
library(broom)

# ---- Load PLINK data ----
ex.data   = read.pedfile(file = "extra.ped", snps = "extra.map")
phenotype = ex.data$fam$affected - 1   # recode 1/2 -> 0/1
genotypes = ex.data$genotypes
snp.ids   = as.character(ex.data$map$snp.names)

cat("Subjects:", nrow(ex.data$fam), "\n")
cat("SNPs:", length(snp.ids), "\n")
cat("Phenotype distribution:\n")
print(table(phenotype))

# Convert genotypes to character data frame
genotypes.df = data.frame(as(genotypes, "character"))
colnames(genotypes.df) = snp.ids

# Example contingency table for one SNP
cat("\nContingency table for rs630969:\n")
print(table(phenotype, genotypes.df$rs630969,
            dnn = c("phenotype", "genotype")))

# ---- Fisher's Exact Test across all SNPs ----
observed.tables.list = sapply(genotypes.df, function(x)
  table(phenotype, x, dnn = c("phenotype", "genotype")))

fish_fn = function(i) {
  cbind(snp.ids[i], fisher.test(observed.tables.list[[i]])$p.value)
}

fish.df = data.frame(t(sapply(1:ncol(genotypes.df), fish_fn)))
colnames(fish.df) = c("rs", "p_value")

fish.results = fish.df %>%
  mutate_at("p_value", as.character) %>%
  mutate_at("p_value", as.numeric) %>%
  arrange(p_value)

cat("\nTop SNPs by Fisher exact test:\n")
print(head(fish.results, 10))

# ---- Logistic regression for top SNP ----
i  = which(snp.ids == fish.results$rs[1])
A1 = ex.data$map$allele.1[i]
A2 = ex.data$map$allele.2[i]
geno.labels = c(paste0(A1, A1), paste0(A1, A2), paste0(A2, A2))

oneSNP.df = data.frame(
  genotypes  = as.factor(genotypes.df[[i]]),
  phenotypes = as.numeric(phenotype)
)

lr.plot = ggplot(oneSNP.df, aes(x = genotypes, y = phenotypes)) +
  geom_point(position = position_jitter(w = 0.1, h = 0.1)) +
  stat_smooth(method = "glm", method.args = list(family = "binomial")) +
  xlim(geno.labels) +
  ggtitle(paste("Logistic regression:", snp.ids[i])) +
  theme_minimal()
print(lr.plot)

# ---- Logistic regression p-values for all SNPs ----
pheno.factor = factor(phenotype, labels = c(0, 1))

lr_all_fn = function(i) {
  lr  = glm(pheno.factor ~ genotypes.df[[i]], family = binomial)
  td  = tidy(lr)
  data.frame(snp = snp.ids[i],
             coef = td$estimate[2],
             pval = td$p.value[2])
}

lr_results = do.call(rbind, lapply(1:ncol(genotypes.df), lr_all_fn))
lr_results = lr_results %>% arrange(pval)

cat("\nTop SNPs by logistic regression:\n")
print(head(lr_results, 10))

# ---- Manhattan plot ----
lr_results$index = 1:nrow(lr_results)
lr_results$logp  = -log10(lr_results$pval)

ggplot(lr_results, aes(x = index, y = logp)) +
  geom_point(color = "steelblue", size = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  xlab("SNP") +
  ylab("-log10(p-value)") +
  ggtitle("Manhattan Plot — SNP-Phenotype Associations") +
  theme_minimal()
