SEA-AD Middle Temporal Gyrus Single-Cell RNA Sequencing Analysis
================================================================

Author: Divya Reddy Konda
Institution: Indiana University
Date: April 2026
GitHub: https://github.com/DivyaReddy262/seaad-mtg-scrna


OVERVIEW
--------

This project performs a comprehensive single-cell RNA sequencing analysis on the
Seattle Alzheimer's Disease Brain Cell Atlas (SEA-AD) Middle Temporal Gyrus dataset.
The goal is to characterize the cellular composition of the human brain and identify
gene expression changes associated with Alzheimer's Disease by comparing Dementia
versus No dementia brain cells at single-cell resolution.


DATASET
-------

Source: AWS S3 - sea-ad-single-cell-profiling
File: SEAAD_MTG_RNAseq_all-nuclei.2024-02-13.h5ad
Size: 43 GB
Format: AnnData (h5ad)
Total nuclei: 1,957,283
Total genes: 36,601
Brain region: Middle Temporal Gyrus (MTG)
Subset used: 50,000 cells for this analysis

Cognitive status in subset:
  Dementia - 20,683 cells (41.4 percent)
  No dementia - 25,806 cells (51.6 percent)
  Reference - 3,511 cells (7.0 percent)


ENVIRONMENT
-----------

HPC: Indiana University Quartz Cluster
Scheduler: SLURM
R version: 4.4.1

Key packages:
  Seurat 5.5.0
  HDF5Array (Bioconductor)
  rhdf5 (Bioconductor)
  clusterProfiler (Bioconductor)
  org.Hs.eg.db (Bioconductor)
  slingshot (Bioconductor)
  speckle (Bioconductor)
  ggplot2
  dplyr


PIPELINE STEPS
--------------

Step 01 - Data Loading

The raw h5ad file was read using H5SparseMatrix from the HDF5Array package, which
handles sparse HDF5 data without loading the full 43 GB into memory. A Seurat object
was created containing the raw UMI count matrix and cell-level metadata including
official SEA-AD cell type labels, donor IDs, and cognitive status.

Step 02 - Quality Control

Filtering thresholds applied:
  nFeature_RNA between 500 and 10,000
  nCount_RNA greater than 500
  frac_mito less than 5 percent

Result: 50,000 input cells reduced to 46,285 cells after filtering.

Step 03 - Normalization and PCA

Log normalization, top 2,000 variable genes, scaling, and PCA (30 PCs).
20 PCs selected based on elbow plot.

Step 04 - UMAP and Clustering

K-nearest neighbor graph using top 20 PCs, Louvain clustering at resolution 0.5,
and UMAP visualization. Result: 27 distinct cell clusters.

Step 05 - Marker Gene Identification

Wilcoxon rank-sum test via FindAllMarkers (min.pct=0.25, logfc.threshold=0.25).

Step 06 - Cell Type Annotation

Official SEA-AD subclass labels from the Allen Institute applied directly,
providing 26 distinct cell subtypes.

Step 07 - Differential Expression Analysis

Comparison: Dementia versus No dementia within each major cell type.
Method: Wilcoxon rank-sum test. Significance: adjusted p-value less than 0.05.

DE genes found per cell type:
  L2/3 IT - 3,025 genes
  Astrocyte - 3,025 genes
  Pvalb - 3,091 genes
  Sst - 3,400 genes
  Microglia-PVM - 2,836 genes
  Vip - 2,297 genes
  Oligodendrocyte - 1,881 genes

Step 08 - Pathway Enrichment Analysis

GO Biological Process and KEGG pathway enrichment using clusterProfiler.
Benjamini-Hochberg FDR correction, q-value cutoff 0.05.

Step 09 - Heatmap

Top 5 DE genes per cell type visualized across all cell populations.
100 cells downsampled per cell type.

Step 10 - Cell Proportion Analysis

Cell type proportions computed for Dementia, No dementia, and Reference groups.

Step 11 - Trajectory Analysis

Slingshot pseudotime on excitatory cortical neuron subtypes.
Root: L2/3 IT neurons. Cell types: L2/3 IT, L4 IT, L5 IT, L5 ET, L5/6 NP, L6 IT, L6 CT, L6b.


KEY FINDINGS
------------

1. L2/3 IT neurons are selectively vulnerable in AD. Highest number of DE genes
   (3,025) and reduced proportion in Dementia donors (18.4 percent vs 22.0 percent).

2. Mitochondrial dysfunction is widespread. MT-ND3 and MTRNR2L12 upregulated
   across multiple neuronal cell types, consistent with the mitochondrial cascade
   hypothesis of Alzheimer's Disease.

3. Astrocytes show reactive astrogliosis. HSPB1 and CHI3L1 upregulated,
   reflecting neuroinflammatory stress response in AD.

4. Microglia show neuroinflammation signatures. Phagocytosis and cytokine
   production pathways enriched, consistent with disease-associated microglial
   phenotype in AD.

5. Oligodendrocytes show altered myelination. PLP1 and myelin-related genes
   differentially expressed, suggesting white matter changes in AD.

6. Cell type proportions largely preserved. AD primarily affects gene expression
   rather than causing overt cell death at this stage.

7. Cortical neuron trajectory runs from L2/3 to L5 to L4. Slingshot pseudotime
   reveals transcriptional hierarchy among excitatory cortical neurons.


REPOSITORY STRUCTURE
--------------------

seaad-mtg-scrna/
  README.md
  scripts/
    load_data.R
    02_normalize_pca.R / .sh
    03_umap_cluster.R / .sh
    04_markers.R / .sh
    05_annotate.R / .sh
    06_de_analysis.R / .sh
    07_pathway_enrichment.R / .sh
    08_heatmap.R / .sh
    09_cell_proportion.R / .sh
    10_trajectory.R / .sh
  results/
    qc/                         QC violin and scatter plots
    clustering/                 Elbow plot, marker gene CSVs
    umap/                       UMAP plots (clusters, celltypes, cogn, annotated)
    differential_expression/
      plots/                    Volcano plots per cell type
      tables/                   DE results and top DE genes CSVs
    pathway_enrichment/
      plots/                    GO and KEGG dot plots per cell type
      tables/                   GO and KEGG result CSVs
    heatmap/                    Heatmap of top DE genes
    cell_proportion/            Barplot, dotplot, proportions CSV
    trajectory/                 Pseudotime and cell type trajectory plots


HOW TO RUN
----------

Module Setup

module load gnu/9.3.0
module load r/4.4.1
export LD_PRELOAD=/N/soft/rhel8/gcc/14.2.0/lib64/libstdc++.so.6
export LD_LIBRARY_PATH=/N/soft/rhel8/gcc/14.2.0/lib64:$LD_LIBRARY_PATH

Download Data

wget https://sea-ad-single-cell-profiling.s3.amazonaws.com/MTG/RNAseq/SEAAD_MTG_RNAseq_all-nuclei.2024-02-13.h5ad

Package Installation (run once inside R)

install.packages("pak")
pak::pkg_install("Seurat")
BiocManager::install("HDF5Array")
BiocManager::install("rhdf5")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("slingshot")
BiocManager::install("SingleCellExperiment")
BiocManager::install("speckle")

Submit Jobs in Order

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

Resource Requirements

Data loading: 128 GB RAM, 8 cores, 6 hours
Normalization and PCA: 256 GB RAM, 8 cores, 12 hours
UMAP and clustering: 128 GB RAM, 8 cores, 6 hours
Differential expression: 128 GB RAM, 8 cores, 8 hours
Pathway enrichment: 64 GB RAM, 4 cores, 4 hours
All other steps: 64 GB RAM, 4 cores, 2 hours

Important Notes

Use H5SparseMatrix instead of zellkonverter to read h5ad files on this cluster.
Always set LD_PRELOAD before starting R to fix GCC library mismatch.
Use -A workshop in SLURM scripts for better queue priority.
Scratch storage has 60-day deletion policy - files are backed up to Slate.


REFERENCES
----------

Leng K, et al. (2021). Molecular characterization of selectively vulnerable neurons
in Alzheimer's disease. Nature Neuroscience, 24(2), 276-287.

Allen Institute for Brain Science. (2024). Seattle Alzheimer's Disease Brain Cell
Atlas. https://sea-ad-single-cell-profiling.s3.amazonaws.com

Hao Y, et al. (2024). Dictionary learning for integrative, multimodal and scalable
single-cell analysis. Nature Biotechnology, 42, 293-304.

Street K, et al. (2018). Slingshot: cell lineage and pseudotime inference for
single-cell transcriptomics. BMC Genomics, 19, 477.

Wu T, et al. (2021). clusterProfiler 4.0: A universal enrichment tool for
interpreting omics data. The Innovation, 2(3), 100141.
