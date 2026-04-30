.libPaths("/N/scratch/dikonda/R_libs")
library(Seurat)
library(ggplot2)
library(dplyr)

cat("Loading clustered object...\n")
seurat_obj <- readRDS("/N/scratch/dikonda/seurat_project/results/seurat_clustered.rds")

cat("Finding all markers...\n")
markers <- FindAllMarkers(seurat_obj,
                          only.pos = TRUE,
                          min.pct = 0.25,
                          logfc.threshold = 0.25)

cat("Top markers per cluster:\n")
top_markers <- markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)
print(top_markers)

write.csv(markers, 
          "/N/scratch/dikonda/seurat_project/results/all_markers.csv",
          row.names = FALSE)
write.csv(top_markers,
          "/N/scratch/dikonda/seurat_project/results/top5_markers.csv",
          row.names = FALSE)
cat("Markers saved!\n")
