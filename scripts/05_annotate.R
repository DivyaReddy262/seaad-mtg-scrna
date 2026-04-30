# ============================================
# Script: 05_annotate.R
# Description: Cell type annotation using
#              official SEA-AD subclass labels
# ============================================

.libPaths("/N/scratch/dikonda/R_libs")
library(Seurat)
library(ggplot2)

cat("Loading clustered object...\n")
seurat_obj <- readRDS("/N/scratch/dikonda/seurat_project/results/seurat_clustered.rds")

# Use official SEA-AD subclass annotations directly
cat("Applying official SEA-AD subclass labels...\n")
Idents(seurat_obj) <- seurat_obj@meta.data$subclass
seurat_obj@meta.data$cell_annotation <- seurat_obj@meta.data$subclass

# Plot with official labels
cat("Saving annotated UMAP...\n")
p <- DimPlot(seurat_obj, reduction = "umap",
             label = TRUE, repel = TRUE,
             label.size = 3, pt.size = 0.3) +
     theme_bw() +
     ggtitle("Official SEA-AD Cell Type Annotations") +
     theme(plot.background = element_rect(fill = "white"),
           panel.background = element_rect(fill = "white")) +
     guides(color = guide_legend(override.aes = list(size = 4)))

ggsave("/N/scratch/dikonda/seurat_project/results/umap_annotated.png",
       p, width = 16, height = 10, bg = "white")

cat("Cell type distribution:\n")
print(table(seurat_obj@meta.data$cell_annotation))

saveRDS(seurat_obj,
        "/N/scratch/dikonda/seurat_project/results/seurat_annotated.rds",
        compress = FALSE)
cat("Done!\n")
