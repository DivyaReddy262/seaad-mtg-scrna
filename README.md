# SEA-AD Middle Temporal Gyrus: Single-Cell RNA Sequencing Analysis

**Author:** Divya Reddy Konda  
**Institution:** Indiana University  
**Date:** April 2026  
**GitHub:** https://github.com/DivyaReddy262/seaad-mtg-scrna

---

## Overview

This project performs a comprehensive single-cell RNA sequencing (scRNA-seq) analysis on the Seattle Alzheimer's Disease Brain Cell Atlas (SEA-AD) Middle Temporal Gyrus (MTG) dataset. The goal is to characterize the cellular composition of the human brain and identify gene expression changes associated with Alzheimer's Disease by comparing Dementia versus No dementia brain cells at single-cell resolution.

---

## Dataset

| Property | Details |
|---|---|
| Source | AWS S3 - sea-ad-single-cell-profiling |
| File | SEAAD_MTG_RNAseq_all-nuclei.2024-02-13.h5ad |
| Size | 43 GB |
| Format | AnnData (h5ad) |
| Total nuclei | 1,957,283 |
| Total genes | 36,601 |
| Brain region | Middle Temporal Gyrus (MTG) |
| Subset used | 50,000 cells |

**Cognitive status breakdown:**

| Group | Cells | Percentage |
|---|---|---|
| Dementia | 20,683 | 41.4% |
| No dementia | 25,806 | 51.6% |
| Reference | 3,511 | 7.0% |

---

## Environment

- **HPC:** Indiana University Quartz Cluster
- **Scheduler:** SLURM
- **R version:** 4.4.1
- **Key packages:** Seurat 5.5.0, HDF5Array, clusterProfiler, slingshot, speckle, ggplot2, dplyr

---

## Pipeline Steps

### Step 01 - Data Loading
The raw h5ad file was read using H5SparseMatrix from the HDF5Array package, which handles sparse HDF5 data without loading the full 43 GB into memory. A Seurat object was created containing the raw UMI count matrix and cell-level metadata including official SEA-AD cell type labels, donor IDs, and cognitive status.

### Step 02 - Quality Control
Three metrics were computed per cell and filtered:
- **nFeature_RNA:** between 500 and 10,000
- **nCount_RNA:** greater than 500
- **frac_mito:** less than 5%

Result: 50,000 cells reduced to 46,285 cells (7.4% removed).

### Step 03 - Normalization and PCA
Log normalization, top 2,000 variable genes, scaling, and PCA (30 PCs). 20 PCs selected based on elbow plot.

### Step 04 - UMAP and Clustering
K-nearest neighbor graph using top 20 PCs, Louvain clustering at resolution 0.5, and UMAP visualization. Result: 27 distinct cell clusters.

### Step 05 - Marker Gene Identification
Wilcoxon rank-sum test via FindAllMarkers (min.pct=0.25, logfc.threshold=0.25).

### Step 06 - Cell Type Annotation
Official SEA-AD subclass labels from the Allen Institute applied directly, providing 26 distinct cell subtypes.

### Step 07 - Differential Expression Analysis
Comparison: Dementia versus No dementia within each major cell type using Wilcoxon rank-sum test (adjusted p-value < 0.05).

| Cell Type | DE Genes |
|---|---|
| L2/3 IT | 3,025 |
| Astrocyte | 3,025 |
| Pvalb | 3,091 |
| Sst | 3,400 |
| Microglia-PVM | 2,836 |
| Vip | 2,297 |
| Oligodendrocyte | 1,881 |

### Step 08 - Pathway Enrichment Analysis
GO Biological Process and KEGG pathway enrichment using clusterProfiler with Benjamini-Hochberg FDR correction (q-value < 0.05).

### Step 09 - Heatmap
Top 5 DE genes per cell type visualized across all cell populations with 100 cells downsampled per cell type.

### Step 10 - Cell Proportion Analysis
Cell type proportions computed for Dementia, No dementia, and Reference groups using speckle.

### Step 11 - Trajectory Analysis
Slingshot pseudotime on excitatory cortical neuron subtypes (L2/3 IT, L4 IT, L5 IT, L5 ET, L5/6 NP, L6 IT, L6 CT, L6b) rooted at L2/3 IT neurons.

---

## Key Findings

1. **L2/3 IT neurons are selectively vulnerable in AD.** Highest number of DE genes (3,025) and reduced proportion in Dementia donors (18.4% vs 22.0%).

2. **Mitochondrial dysfunction is widespread.** MT-ND3 and MTRNR2L12 upregulated across multiple neuronal cell types, consistent with the mitochondrial cascade hypothesis.

3. **Astrocytes show reactive astrogliosis.** HSPB1 and CHI3L1 upregulated, reflecting neuroinflammatory stress response in AD.

4. **Microglia show neuroinflammation signatures.** Phagocytosis and cytokine production pathways enriched, consistent with disease-associated microglial phenotype.

5. **Oligodendrocytes show altered myelination.** PLP1 and myelin-related genes differentially expressed, suggesting white matter changes in AD.

6. **Cell type proportions largely preserved.** AD primarily affects gene expression rather than causing overt cell death at this stage.

7. **Cortical neuron trajectory runs from L2/3 to L5 to L4.** Slingshot pseudotime reveals transcriptional hierarchy among excitatory cortical neurons.

---

## Repository Structure
## Repository Structure
seaad-mtg-scrna/
├── README.md
├── scripts/
│   ├── load_data.R
│   ├── 02_normalize_pca.R / .sh
│   ├── 03_umap_cluster.R / .sh
│   ├── 04_markers.R / .sh
│   ├── 05_annotate.R / .sh
│   ├── 06_de_analysis.R / .sh
│   ├── 07_pathway_enrichment.R / .sh
│   ├── 08_heatmap.R / .sh
│   ├── 09_cell_proportion.R / .sh
│   └── 10_trajectory.R / .sh
└── results/
├── qc/                        QC violin and scatter plots
├── clustering/                Elbow plot, marker gene CSVs
├── umap/                      UMAP plots
├── differential_expression/
│   ├── plots/                 Volcano plots per cell type
│   └── tables/                DE results CSVs
├── pathway_enrichment/
│   ├── plots/                 GO and KEGG dot plots
│   └── tables/                GO and KEGG result CSVs
├── heatmap/                   Heatmap of top DE genes
├── cell_proportion/           Proportion plots and CSV
└── trajectory/                Pseudotime plots
---

## How to Run

**Module Setup**
```bash
module load gnu/9.3.0
module load r/4.4.1
export LD_PRELOAD=/N/soft/rhel8/gcc/14.2.0/lib64/libstdc++.so.6
export LD_LIBRARY_PATH=/N/soft/rhel8/gcc/14.2.0/lib64:$LD_LIBRARY_PATH
```

**Download Data**
```bash
wget https://sea-ad-single-cell-profiling.s3.amazonaws.com/MTG/RNAseq/SEAAD_MTG_RNAseq_all-nuclei.2024-02-13.h5ad
```

**Package Installation**
```r
install.packages("pak")
pak::pkg_install("Seurat")
BiocManager::install(c("HDF5Array", "rhdf5", "clusterProfiler",
                       "org.Hs.eg.db", "slingshot", "SingleCellExperiment", "speckle"))
```

**Submit Jobs in Order**
```bash
sbatch scripts/load_data.sh
sbatch scripts/02_normalize_pca.sh
sbatch scripts/03_umap_cluster.sh
sbatch scripts/04_markers.sh
sbatch scripts/05_annotate.sh
sbatch scripts/06_de_analysis.sh
sbatch scripts/07_pathway_enrichment.sh
sbatch scripts/08_heatmap.sh
sbatch scripts/09_cell_proportion.sh
sbatch scripts/10_trajectory.sh
```

**Resource Requirements**

| Step | RAM | Cores | Time |
|---|---|---|---|
| Data loading | 128 GB | 8 | 6 hrs |
| Normalization and PCA | 256 GB | 8 | 12 hrs |
| UMAP and clustering | 128 GB | 8 | 6 hrs |
| Differential expression | 128 GB | 8 | 8 hrs |
| Pathway enrichment | 64 GB | 4 | 4 hrs |
| All other steps | 64 GB | 4 | 2 hrs |

**Important Notes**
- Use H5SparseMatrix instead of zellkonverter to read h5ad files on this cluster
- Always set LD_PRELOAD before starting R to fix GCC library mismatch
- Use `-A workshop` in SLURM scripts for better queue priority
- Scratch storage has 60-day deletion policy — files are backed up to Slate

---

## References

Leng K, et al. (2021). Molecular characterization of selectively vulnerable neurons in Alzheimer's disease. *Nature Neuroscience*, 24(2), 276-287.

Allen Institute for Brain Science. (2024). Seattle Alzheimer's Disease Brain Cell Atlas. https://sea-ad-single-cell-profiling.s3.amazonaws.com

Hao Y, et al. (2024). Dictionary learning for integrative, multimodal and scalable single-cell analysis. *Nature Biotechnology*, 42, 293-304.

Street K, et al. (2018). Slingshot: cell lineage and pseudotime inference for single-cell transcriptomics. *BMC Genomics*, 19, 477.

Wu T, et al. (2021). clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. *The Innovation*, 2(3), 100141.
