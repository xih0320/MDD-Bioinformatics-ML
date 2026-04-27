# ===============================
# 05_STRING_PPI_network.R
# STRING PPI Network Analysis
# ===============================

# 1. Load required packages
library(STRINGdb)
library(dplyr)
library(ggplot2)
library(igraph)

# 2. Create output folders
dir.create("results", showWarnings = FALSE)
dir.create("plots", showWarnings = FALSE)

# 3. Define file paths
sig_gene_file <- "results/MDD_DE_genes_sig.csv"

mapped_string_file <- "results/STRING_mapped_genes.csv"
ppi_edge_file <- "results/STRING_PPI_edges.csv"
hub_gene_file <- "results/STRING_hub_genes.csv"

ppi_plot_file <- "plots/STRING_PPI_network.png"
hub_plot_file <- "plots/STRING_hub_genes_barplot.png"

# 4. Load significant DE genes
de_genes <- read.csv(sig_gene_file, check.names = FALSE)

print(colnames(de_genes))
head(de_genes)

# 5. Extract gene symbols
gene_symbols <- de_genes$`Gene Symbol`

gene_symbols <- gene_symbols[!is.na(gene_symbols)]
gene_symbols <- gene_symbols[gene_symbols != ""]
gene_symbols <- unique(gene_symbols)

cat("Number of unique gene symbols:", length(gene_symbols), "\n")

gene_df <- data.frame(
  gene = gene_symbols,
  stringsAsFactors = FALSE
)

# 6. Initialize STRING database
string_db <- STRINGdb$new(
  version = "12.0",
  species = 9606,
  score_threshold = 400,
  input_directory = ""
)

# 7. Map gene symbols to STRING IDs
mapped_genes <- string_db$map(
  gene_df,
  "gene",
  removeUnmappedRows = TRUE
)

write.csv(
  mapped_genes,
  mapped_string_file,
  row.names = FALSE
)

cat("Number of genes mapped to STRING:", nrow(mapped_genes), "\n")

# 8. Retrieve PPI interactions
ppi_edges <- string_db$get_interactions(mapped_genes$STRING_id)

write.csv(
  ppi_edges,
  ppi_edge_file,
  row.names = FALSE
)

cat("Number of PPI edges:", nrow(ppi_edges), "\n")

# 9. Build igraph network
ppi_graph <- graph_from_data_frame(
  ppi_edges[, c("from", "to")],
  directed = FALSE
)

# 10. Calculate node degree
degree_values <- degree(ppi_graph)

hub_genes <- data.frame(
  STRING_id = names(degree_values),
  Degree = as.numeric(degree_values),
  stringsAsFactors = FALSE
)

# 11. Add gene symbols back to hub table
hub_genes <- hub_genes %>%
  left_join(
    mapped_genes %>% select(STRING_id, gene),
    by = "STRING_id"
  ) %>%
  arrange(desc(Degree))

write.csv(
  hub_genes,
  hub_gene_file,
  row.names = FALSE
)

cat("Top hub genes:\n")
print(head(hub_genes, 10))

# 12. Plot PPI network
png(ppi_plot_file, width = 1000, height = 800)

plot(
  ppi_graph,
  vertex.size = 8,
  vertex.label = hub_genes$gene[match(V(ppi_graph)$name, hub_genes$STRING_id)],
  vertex.label.cex = 0.8,
  edge.width = 1,
  main = "STRING PPI Network of Significant DE Genes"
)

dev.off()

cat("PPI network plot saved to:", ppi_plot_file, "\n")

# 13. Plot top hub genes
top_hub_genes <- hub_genes %>%
  filter(!is.na(gene)) %>%
  head(10)

png(hub_plot_file, width = 1000, height = 700)

ggplot(top_hub_genes, aes(x = reorder(gene, Degree), y = Degree)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    title = "Top Hub Genes in STRING PPI Network",
    x = "Gene Symbol",
    y = "Degree"
  ) +
  theme_minimal()

dev.off()

cat("Hub gene barplot saved to:", hub_plot_file, "\n")

cat("\nSTRING PPI network analysis completed successfully.\n")
