.libPaths("/N/scratch/dikonda/R_libs")
library(Seurat)
library(ggplot2)
library(dplyr)

cat("Loading annotated object...\n")
seurat_obj <- readRDS("/N/scratch/dikonda/seurat_project/results/seurat_annotated.rds")

# Set identity to subclass for cell-type specific DE
Idents(seurat_obj) <- seurat_obj@meta.data$subclass

cat("Cell type distribution:\n")
print(table(seurat_obj@meta.data$subclass))
print(table(seurat_obj@meta.data$cogn_status))

# Focus on major cell types with enough cells
# Compare Dementia vs No dementia within each cell type
cell_types <- c("L2/3 IT", "L5 IT", "L4 IT", "Astrocyte", 
                "Oligodendrocyte", "Microglia-PVM", "Vip", 
                "Pvalb", "Sst")

de_results <- list()

for (ct in cell_types) {
  cat("\nRunning DE for:", ct, "\n")
  
  # Subset to this cell type
  cells_ct <- subset(seurat_obj, subclass == ct)
  
  # Check enough cells in each group
  n_dementia <- sum(cells_ct@meta.data$cogn_status == "Dementia")
  n_normal   <- sum(cells_ct@meta.data$cogn_status == "No dementia")
  
  cat("  Dementia cells:", n_dementia, "\n")
  cat("  No dementia cells:", n_normal, "\n")
  
  if (n_dementia < 10 | n_normal < 10) {
    cat("  Skipping - not enough cells\n")
    next
  }
  
  # Set identity to cognitive status
  Idents(cells_ct) <- cells_ct@meta.data$cogn_status
  
  # Run DE: Dementia vs No dementia
  tryCatch({
    de <- FindMarkers(cells_ct,
                      ident.1 = "Dementia",
                      ident.2 = "No dementia",
                      min.pct = 0.1,
                      logfc.threshold = 0.25)
    de$gene <- rownames(de)
    de$cell_type <- ct
    de_results[[ct]] <- de
    cat("  Found", nrow(de), "DE genes\n")
  }, error = function(e) {
    cat("  Error:", e$message, "\n")
  })
}

# Combine all results
cat("\nCombining results...\n")
all_de <- bind_rows(de_results)

# Save results
write.csv(all_de,
          "/N/scratch/dikonda/seurat_project/results/de_results.csv",
          row.names = FALSE)

# Top DE genes per cell type
top_de <- all_de %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cell_type) %>%
  top_n(n = 10, wt = abs(avg_log2FC))

write.csv(top_de,
          "/N/scratch/dikonda/seurat_project/results/top_de_genes.csv",
          row.names = FALSE)

cat("\nTop DE genes per cell type:\n")
print(top_de %>% select(cell_type, gene, avg_log2FC, p_val_adj))

# Plot top DE genes as a heatmap
cat("\nSaving volcano plots...\n")
for (ct in names(de_results)) {
  de <- de_results[[ct]]
  de$gene <- rownames(de)
  de$significant <- de$p_val_adj < 0.05 & abs(de$avg_log2FC) > 0.25
  de$direction <- ifelse(de$avg_log2FC > 0, "Up in AD", "Down in AD")

  p <- ggplot(de, aes(x = avg_log2FC, y = -log10(p_val_adj),
                       color = significant)) +
       geom_point(alpha = 0.6, size = 1) +
       scale_color_manual(values = c("grey", "red")) +
       geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed") +
       geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
       theme_bw() +
       ggtitle(paste("DE:", ct, "- Dementia vs No dementia")) +
       xlab("log2 Fold Change") +
       ylab("-log10 adjusted p-value") +
       theme(plot.background = element_rect(fill = "white"),
             panel.background = element_rect(fill = "white"))

  filename <- paste0("/N/scratch/dikonda/seurat_project/results/volcano_",
                     gsub("/", "_", ct), ".png")
  ggsave(filename, p, width = 8, height = 6, bg = "white")
}

cat("Done! DE analysis complete.\n")
