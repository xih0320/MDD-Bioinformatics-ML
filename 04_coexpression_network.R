# ============================================================
# MDD Gene Expression Analysis
# 04: Gene Co-expression Network Analysis
# ============================================================
# Requires: GxS_covfilter from 01_data_preprocessing.R
# Uses stricter CoV filter for manageable network size

library(preprocessCore)
library(dplyr)
library(RMThreshold)
library(igraph)

# Stricter CoV filter for network analysis
cov_thresh = .05
CoV_values = apply(mddExprData_quantileLog2, 1, function(x) { sd(x) / abs(mean(x)) })
sd_values  = apply(mddExprData_quantileLog2, 1, function(x) { sd(x) })

GxS_covfilter = mddExprData_quantileLog2[CoV_values < cov_thresh & sd_values > 0, ]
cat("Genes for network analysis:", nrow(GxS_covfilter), "\n")

# ---- Build co-expression matrix ----
mddCor = cor(t(GxS_covfilter))

# ---- Random Matrix Theory threshold to remove noise edges ----
res     = rm.get.threshold(mddCor, interval = c(0.4, 0.7), nr.thresholds = 10)
cleaned = rm.denoise.mat(mddCor, 0.5, keep.diag = FALSE)

# Binary adjacency matrix
adjMat = cleaned != 0
rownames(adjMat) = rownames(GxS_covfilter)
colnames(adjMat) = rownames(GxS_covfilter)

# ---- Build graph and remove isolated nodes ----
ig        = graph.adjacency(adjMat, mode = "undirected")
igDegrees = rowSums(adjMat)

hist(igDegrees, col = "steelblue", main = "Degree Distribution (all nodes)",
     xlab = "Degree", ylab = "Frequency")

# Remove degree-0 nodes
adjMatConnected = adjMat[igDegrees != 0, igDegrees != 0]
igConnected     = graph.adjacency(adjMatConnected, mode = "undirected")

igDegreesConnected = rowSums(adjMatConnected)
cat("Connected nodes:", nrow(adjMatConnected), "\n")

par(mfrow = c(1, 2))
plot(igConnected, vertex.size = 1, vertex.label = NA, edge.width = 1,
     main = "Co-expression Network")
hist(igDegreesConnected, col = "steelblue",
     main = "Degree Distribution (connected)", xlab = "Degree")
par(mfrow = c(1, 1))

# ---- Community detection ----

# 1. Fast-greedy (greedy modularity)
igc_clusts     = fastgreedy.community(igConnected)
igc_membership = membership(igc_clusts)
cat("\nFast-greedy clusters:", length(igc_clusts), "\n")
print(sizes(igc_clusts)[1:5])

# 2. Louvain
louvain_clusts = cluster_louvain(igConnected, resolution = 1)
cat("Louvain clusters:", length(louvain_clusts), "\n")

# 3. Leiden
r = quantile(strength(igConnected))[2] / (gorder(igConnected) - 1)
leiden_clusts = cluster_leiden(igConnected, resolution_parameter = r)
cat("Leiden clusters:", length(leiden_clusts), "\n")

# Color network by fast-greedy membership
color_palette = rainbow(length(igc_clusts))
igc_colors    = color_palette[igc_membership]

plot(igConnected, vertex.size = 2, vertex.label = NA,
     vertex.color = igc_colors, edge.width = 0.5,
     main = "Co-expression Network — Fast-Greedy Modules")

# ---- Export top two cluster gene lists for pathway enrichment ----
clust1_genes = names(igc_membership)[igc_membership == 1]
clust2_genes = names(igc_membership)[igc_membership == 2]

write.table(clust1_genes, file = "clust1.txt", row.names = FALSE,
            col.names = FALSE, quote = FALSE)
write.table(clust2_genes, file = "clust2.txt", row.names = FALSE,
            col.names = FALSE, quote = FALSE)

cat("\nCluster 1 size:", length(clust1_genes), "\n")
cat("Cluster 2 size:", length(clust2_genes), "\n")
cat("Gene lists saved to clust1.txt and clust2.txt\n")
