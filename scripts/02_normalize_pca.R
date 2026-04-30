.libPaths("/N/scratch/dikonda/R_libs")
library(Seurat)
library(ggplot2)

cat("Loading filtered object...\n")
seurat_obj <- readRDS("/N/scratch/dikonda/seurat_project/results/seurat_filtered.rds")

cat("Normalizing...\n")
seurat_obj <- NormalizeData(seurat_obj)

cat("Finding variable features...\n")
seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000)

cat("Scaling...\n")
seurat_obj <- ScaleData(seurat_obj)

cat("Running PCA...\n")
seurat_obj <- RunPCA(seurat_obj)

cat("Saving elbow plot...\n")
p <- ElbowPlot(seurat_obj, ndims = 30) +
     theme_bw() +
     theme(panel.background = element_rect(fill = "white"),
           plot.background = element_rect(fill = "white"))

ggsave("/N/scratch/dikonda/seurat_project/results/elbow_plot.png",
       p, width = 8, height = 6, bg = "white")

cat("Saving object...\n")
saveRDS(seurat_obj,
        "/N/scratch/dikonda/seurat_project/results/seurat_pca.rds",
        compress = FALSE)
cat("Done!\n")
