.libPaths("/N/scratch/dikonda/R_libs")
library(HDF5Array)
library(rhdf5)
library(Matrix)
library(Seurat)

h5file <- "/N/scratch/dikonda/seurat_project/data/SEAAD_MTG_RNAseq_all-nuclei.2024-02-13.h5ad"

cat("Loading gene and cell names...\n")
genes <- h5read(h5file, "/var/_index")
cells <- h5read(h5file, "/obs/bc")
cells_unique <- make.unique(as.character(cells))

cat("Loading count matrix...\n")
counts <- H5SparseMatrix(h5file, "layers/UMIs")

cat("Subsetting to 50k cells...\n")
counts_sub <- counts[, 1:50000]
counts_sparse <- as(counts_sub, "CsparseMatrix")
rownames(counts_sparse) <- genes
colnames(counts_sparse) <- cells_unique[1:50000]

cat("Reading category labels...\n")
cell_type_cats <- h5read(h5file, "/obs/__categories/Class")
subclass_cats  <- h5read(h5file, "/obs/__categories/Subclass")
supertype_cats <- h5read(h5file, "/obs/__categories/Supertype")
cogn_cats      <- h5read(h5file, "/obs/__categories/Cognitive Status")
sex_cats       <- h5read(h5file, "/obs/__categories/Sex")

cat("Building metadata...\n")
meta <- data.frame(
  row.names    = cells_unique[1:50000],
  cell_type    = cell_type_cats[as.numeric(h5read(h5file, "/obs/Class")[1:50000]) + 1],
  subclass     = subclass_cats[as.numeric(h5read(h5file, "/obs/Subclass")[1:50000]) + 1],
  supertype    = supertype_cats[as.numeric(h5read(h5file, "/obs/Supertype")[1:50000]) + 1],
  donor_id     = as.character(h5read(h5file, "/obs/Donor ID")[1:50000]),
  cogn_status  = cogn_cats[as.numeric(h5read(h5file, "/obs/Cognitive Status")[1:50000]) + 1],
  sex          = sex_cats[as.numeric(h5read(h5file, "/obs/Sex")[1:50000]) + 1],
  age_at_death = h5read(h5file, "/obs/Age at Death")[1:50000],
  braak        = h5read(h5file, "/obs/Braak")[1:50000],
  n_umis       = h5read(h5file, "/obs/Number of UMIs")[1:50000],
  n_genes      = h5read(h5file, "/obs/Genes detected")[1:50000],
  frac_mito    = h5read(h5file, "/obs/Fraction mitochondrial UMIs")[1:50000]
)

cat("Creating Seurat object...\n")
seurat_obj <- CreateSeuratObject(
  counts       = counts_sparse,
  meta.data    = meta,
  project      = "SEAAD_MTG",
  min.cells    = 3,
  min.features = 200
)

cat("Done!\n")
print(seurat_obj)
cat("\nCell type distribution:\n")
print(table(seurat_obj@meta.data$cell_type))
cat("\nCognitive status:\n")
print(table(seurat_obj@meta.data$cogn_status))

saveRDS(seurat_obj,
        "/N/scratch/dikonda/seurat_project/results/seurat_raw.rds",
        compress = FALSE)
cat("Saved!\n")
