.libPaths("/N/scratch/dikonda/R_libs")
library(Seurat)
library(ggplot2)

cat("Loading PCA object...\n")
seurat_obj <- readRDS("/N/scratch/dikonda/seurat_project/results/seurat_pca.rds")

cat("Finding neighbors...\n")
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)

cat("Finding clusters...\n")
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

cat("Running UMAP...\n")
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)

cat("Saving UMAP plots...\n")
# UMAP by cluster
p1 <- DimPlot(seurat_obj, reduction = "umap", label = TRUE) +
      theme_bw() + ggtitle("Clusters")
ggsave("/N/scratch/dikonda/seurat_project/results/umap_clusters.png",
       p1, width = 8, height = 7, bg = "white")

# UMAP by cell type
p2 <- DimPlot(seurat_obj, reduction = "umap",
              group.by = "cell_type", label = TRUE) +
      theme_bw() + ggtitle("Cell Types")
ggsave("/N/scratch/dikonda/seurat_project/results/umap_celltypes.png",
       p2, width = 10, height = 7, bg = "white")

# UMAP by cognitive status
p3 <- DimPlot(seurat_obj, reduction = "umap",
              group.by = "cogn_status") +
      theme_bw() + ggtitle("Cognitive Status")
ggsave("/N/scratch/dikonda/seurat_project/results/umap_cogn.png",
       p3, width = 8, height = 7, bg = "white")

cat("Saving object...\n")
saveRDS(seurat_obj,
        "/N/scratch/dikonda/seurat_project/results/seurat_clustered.rds",
        compress = FALSE)
cat("Done!\n")
print(seurat_obj)
