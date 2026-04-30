.libPaths("/N/scratch/dikonda/R_libs")
library(Seurat)
library(ggplot2)

cat("Loading clustered object...\n")
seurat_obj <- readRDS("/N/scratch/dikonda/seurat_project/results/seurat_clustered.rds")

# Annotate clusters
new_labels <- c(
  "0"  = "Glutamatergic (L2-3)",
  "1"  = "Oligodendrocyte",
  "2"  = "Excitatory Neuron",
  "3"  = "GABAergic Neuron",
  "4"  = "Fibroblast",
  "5"  = "Glutamatergic Neuron",
  "6"  = "VIP Interneuron",
  "7"  = "Astrocyte",
  "8"  = "PV Interneuron",
  "9"  = "SST Interneuron",
  "10" = "Microglia",
  "11" = "GABAergic Interneuron",
  "12" = "OPC",
  "13" = "Glutamatergic Neuron",
  "14" = "Mast Cell",
  "15" = "Glutamatergic Neuron",
  "16" = "Pericyte",
  "17" = "Glutamatergic Neuron",
  "18" = "GABAergic Neuron",
  "19" = "Endothelial Cell",
  "20" = "Glutamatergic Neuron",
  "21" = "Immune Cell",
  "22" = "Microglia",
  "23" = "Pericyte",
  "24" = "Endothelial Cell",
  "25" = "Fibroblast",
  "26" = "Astrocyte"
)

seurat_obj <- RenameIdents(seurat_obj, new_labels)
seurat_obj@meta.data$cell_annotation <- Idents(seurat_obj)

# Plot annotated UMAP
p <- DimPlot(seurat_obj, reduction = "umap",
             label = TRUE, repel = TRUE,
             label.size = 3) +
     theme_bw() +
     ggtitle("Annotated Cell Types") +
     theme(plot.background = element_rect(fill = "white"),
           panel.background = element_rect(fill = "white"))

ggsave("/N/scratch/dikonda/seurat_project/results/umap_annotated.png",
       p, width = 12, height = 8, bg = "white")

cat("Cell type distribution:\n")
print(table(seurat_obj@meta.data$cell_annotation))

saveRDS(seurat_obj,
        "/N/scratch/dikonda/seurat_project/results/seurat_annotated.rds",
        compress = FALSE)
cat("Done!\n")
