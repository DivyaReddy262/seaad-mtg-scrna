.libPaths("/N/scratch/dikonda/R_libs")
options(expressions = 50000)

library(Seurat)
library(SingleCellExperiment)

cat("Loading annotated object...\n")
seurat_obj <- readRDS("/N/scratch/dikonda/seurat_project/results/seurat_annotated.rds")

# Focus on neuronal cells
cat("Subsetting neuronal cells...\n")
neuronal <- subset(seurat_obj,
                   subclass %in% c("L2/3 IT", "L4 IT", "L5 IT",
                                   "L5/6 NP", "L6 IT", "L6 CT",
                                   "L6b", "L5 ET"))

# Downsample
neuronal <- subset(neuronal,
                   cells = WhichCells(neuronal, downsample = 200))
cat("Cells:", ncol(neuronal), "\n")

# Convert to SCE
sce <- as.SingleCellExperiment(neuronal)

# Load slingshot AFTER Seurat to avoid conflicts
library(slingshot)

# Use full namespace to avoid recursion
cat("Running slingshot...\n")
sce <- slingshot::slingshot(sce,
                            clusterLabels = colData(sce)$subclass,
                            reducedDim    = "UMAP",
                            start.clus    = "L2/3 IT")

pt <- slingshot::slingPseudotime(sce)
cat("Pseudotime calculated!\n")

# Plot
library(ggplot2)
umap_coords <- as.data.frame(reducedDim(sce, "UMAP"))
colnames(umap_coords) <- c("umap_1", "umap_2")
umap_coords$pseudotime <- pt[, 1]
umap_coords$subclass   <- colData(sce)$subclass

p1 <- ggplot(umap_coords, aes(x = umap_1, y = umap_2, color = pseudotime)) +
      geom_point(size = 1, alpha = 0.7) +
      scale_color_gradientn(colors = c("blue", "green", "yellow", "red"),
                            na.value = "grey") +
      theme_bw() +
      ggtitle("Neuronal Trajectory - Pseudotime") +
      theme(plot.background  = element_rect(fill = "white"),
            panel.background = element_rect(fill = "white"))

ggsave("/N/scratch/dikonda/seurat_project/results/trajectory_pseudotime.png",
       p1, width = 10, height = 8, bg = "white")

p2 <- ggplot(umap_coords, aes(x = umap_1, y = umap_2, color = subclass)) +
      geom_point(size = 1, alpha = 0.7) +
      theme_bw() +
      ggtitle("Neuronal Trajectory - Cell Types") +
      theme(plot.background  = element_rect(fill = "white"),
            panel.background = element_rect(fill = "white"))

ggsave("/N/scratch/dikonda/seurat_project/results/trajectory_celltypes.png",
       p2, width = 10, height = 8, bg = "white")

cat("Trajectory analysis done!\n")
