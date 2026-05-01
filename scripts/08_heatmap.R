.libPaths("/N/scratch/dikonda/R_libs")
library(Seurat)
library(ggplot2)
library(dplyr)

cat("Loading annotated object...\n")
seurat_obj <- readRDS("/N/scratch/dikonda/seurat_project/results/seurat_annotated.rds")

# Load top DE genes
top_de <- read.csv("/N/scratch/dikonda/seurat_project/results/de_results.csv") %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cell_type) %>%
  top_n(n = 5, wt = abs(avg_log2FC)) %>%
  ungroup()

top_genes <- unique(top_de$gene)
cat("Total genes for heatmap:", length(top_genes), "\n")

# Set identity to subclass
Idents(seurat_obj) <- seurat_obj@meta.data$subclass

# Downsample for heatmap
seurat_sub <- subset(seurat_obj,
                     subclass %in% unique(top_de$cell_type))
seurat_sub <- subset(seurat_sub,
                     cells = WhichCells(seurat_sub, downsample = 100))

cat("Cells in heatmap:", ncol(seurat_sub), "\n")

seurat_sub <- NormalizeData(seurat_sub)

p <- DoHeatmap(seurat_sub,
               features = top_genes,
               group.by = "subclass",
               size = 3,
               angle = 45) +
     scale_fill_gradientn(colors = c("navy", "white", "red")) +
     theme(axis.text.y = element_text(size = 6),
           plot.background = element_rect(fill = "white"),
           panel.background = element_rect(fill = "white"))

ggsave("/N/scratch/dikonda/seurat_project/results/heatmap_top_de_genes.png",
       p, width = 16, height = 12, bg = "white")
cat("Heatmap saved!\n")
