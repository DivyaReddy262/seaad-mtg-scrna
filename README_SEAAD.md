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

Step 01 - Data Loading (load_data.R)

The raw h5ad file was read using H5SparseMatrix from the HDF5Array package, which
handles sparse HDF5 data without loading the full 43 GB into memory. Gene names,
cell barcodes, and metadata were extracted directly from the HDF5 structure. A
Seurat object was created containing the raw UMI count matrix and cell-level metadata
including official SEA-AD cell type labels, donor IDs, cognitive status, age at death,
sex, Braak stage, and QC metrics.


Step 02 - Quality Control (02_quality_control.Rmd)

Three standard metrics were computed per cell:
  nFeature_RNA: number of unique genes detected per cell
  nCount_RNA: total UMI counts per cell
  frac_mito: fraction of reads mapping to mitochondrial genes

Filtering thresholds applied:
  nFeature_RNA between 500 and 10,000
  nCount_RNA greater than 500
  frac_mito less than 5 percent

Result: 50,000 input cells reduced to 46,285 cells after filtering.
Cells removed: 3,715 (7.4 percent removal rate).


Step 03 - Normalization and PCA (02_normalize_pca.R)

Normalization: Log normalization scaling each cell to 10,000 counts then log-transforming.
Variable features: Top 2,000 most variable genes identified using variance-stabilizing transformation.
Scaling: All variable genes centered and scaled to mean zero and unit variance.
PCA: Top 30 principal components computed on scaled variable genes.
PC selection: 20 PCs selected based on elbow plot where variance explained plateaus.


Step 04 - UMAP and Clustering (03_umap_cluster.R)

Neighbor graph: K-nearest neighbor graph built using top 20 PCs via FindNeighbors.
Clustering: Louvain community detection at resolution 0.5 via FindClusters.
UMAP: 2D visualization computed on top 20 PCs via RunUMAP.
Result: 27 distinct cell clusters identified (numbered 0 through 26).


Step 05 - Marker Gene Identification (04_markers.R)

Method: Wilcoxon rank-sum test via FindAllMarkers.
Parameters: min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE.
Output: Top marker genes per cluster saved to CSV files for annotation.


Step 06 - Cell Type Annotation (05_annotate.R)

Official SEA-AD subclass labels from the Allen Institute were applied directly to
the Seurat object. These expert-curated annotations provide 26 distinct cell subtypes
validated by the SEA-AD consortium, covering cortical layer-specific excitatory neurons,
multiple inhibitory interneuron subtypes, and non-neuronal populations.


Step 07 - Differential Expression Analysis (06_de_analysis.R)

Comparison: Dementia versus No dementia cells within each major cell type.
Method: Wilcoxon rank-sum test via FindMarkers.
Significance threshold: adjusted p-value less than 0.05.

Cell types analyzed and DE genes found:
  L2/3 IT neurons - 3,025 DE genes
  Astrocyte - 3,025 DE genes
  Pvalb interneurons - 3,091 DE genes
  Sst interneurons - 3,400 DE genes
  Microglia-PVM - 2,836 DE genes
  Vip interneurons - 2,297 DE genes
  Oligodendrocyte - 1,881 DE genes


Step 08 - Pathway Enrichment Analysis (07_pathway_enrichment.R)

Method: clusterProfiler with org.Hs.eg.db database.
Tests performed: GO Biological Process enrichment and KEGG pathway enrichment.
Multiple testing correction: Benjamini-Hochberg FDR, q-value cutoff 0.05.
Output: Dot plots and CSV tables saved for each cell type.


Step 09 - Heatmap (08_heatmap.R)

Top 5 DE genes per cell type were visualized in a heatmap across all analyzed cell
populations. Cells were downsampled to 100 per cell type. Expression is shown on
a blue-white-red gradient from low to high expression.


Step 10 - Cell Proportion Analysis (09_cell_proportion.R)

Cell type proportions were computed for Dementia, No dementia, and Reference groups.
Results visualized as a stacked bar plot and a dot-line proportion change plot.


Step 11 - Trajectory Analysis (10_trajectory.R)

Method: Slingshot pseudotime on excitatory cortical neuron subtypes.
Cell types included: L2/3 IT, L4 IT, L5 IT, L5 ET, L5/6 NP, L6 IT, L6 CT, L6b.
Root cell type: L2/3 IT neurons (most vulnerable in AD, superficial cortical layer).
Output: Pseudotime trajectory plots colored by cell type and pseudotime value.


RESULTS SUMMARY
---------------

Cell Types Identified

26 distinct cell subtypes were annotated using official SEA-AD expert labels.

Excitatory neurons (Glutamatergic):
  L2/3 IT - 9,442 cells - superficial cortical output neurons
  L5 IT - 6,581 cells - middle layer intratelencephalic neurons
  L4 IT - 4,779 cells - thalamorecipient layer neurons
  L6 IT - 1,132 cells - deep layer neurons
  L6b - 391 cells - subplate neurons
  L6 CT - 488 cells - corticothalamic projection neurons
  L5 ET - 46 cells - extratelencephalic projection neurons
  L5/6 NP - 776 cells - near-projecting neurons
  L6 IT Car3 - 710 cells

Inhibitory neurons (GABAergic):
  Vip - 3,098 cells - vasoactive intestinal peptide interneurons
  Sst - 1,798 cells - somatostatin interneurons
  Pvalb - 2,388 cells - parvalbumin fast-spiking interneurons
  Lamp5 - 1,341 cells
  Chandelier - 337 cells
  Sncg - 756 cells
  Lamp5 Lhx6 - 709 cells
  Sst Chodl - 37 cells
  Pax6 - 316 cells

Non-neuronal cells:
  Oligodendrocyte - 3,696 cells - myelin-producing glia
  Astrocyte - 2,552 cells - support and homeostasis
  VLMC - 2,009 cells - vascular leptomeningeal cells
  Microglia-PVM - 1,610 cells - brain immune cells
  OPC - 1,208 cells - oligodendrocyte precursors
  Endothelial - 85 cells


Key Scientific Findings

Finding 1: L2/3 IT neurons are selectively vulnerable in Alzheimer's Disease.
These superficial cortical neurons showed the highest number of DE genes at 3,025,
including upregulation of mitochondrial dysfunction markers MT-ND3 and MTRNR2L12.
Their proportion in Dementia donors (18.4 percent) was notably lower than in No
dementia donors (22.0 percent), consistent with selective neuronal loss in early AD.

Finding 2: Mitochondrial dysfunction is a widespread AD signature.
Mitochondrial genes MT-ND3 and MTRNR2L12 were upregulated across multiple neuronal
cell types in AD, consistent with the mitochondrial cascade hypothesis of Alzheimer's
Disease. This suggests energy metabolism failure as an early pathological event.

Finding 3: Astrocytes show reactive astrogliosis in AD.
Upregulation of HSPB1 (heat shock protein) and CHI3L1 (neuroinflammation marker)
in astrocytes reflects activation of stress response and neuroinflammatory pathways,
consistent with reactive astrogliosis, a hallmark of neurodegeneration.

Finding 4: Microglia show neuroinflammation signatures.
Enrichment of phagocytosis, cytokine production, and immune response pathways in
Microglia-PVM cells is consistent with the disease-associated microglial phenotype
observed in Alzheimer's Disease.

Finding 5: Oligodendrocytes show altered myelination gene expression.
PLP1 and other myelin-related genes show differential expression in AD
oligodendrocytes, suggesting white matter changes contribute to cognitive dysfunction.

Finding 6: Cell type proportions are largely preserved.
The cell proportion analysis shows broadly similar cellular composition across
cognitive status groups, indicating that AD primarily affects gene expression rather
than causing overt cell death at this stage in the disease course.

Finding 7: Cortical neuron trajectory runs from L2/3 to L5 to L4.
Slingshot pseudotime analysis reveals a transcriptional hierarchy among excitatory
cortical neurons anchored at L2/3 IT neurons, providing a framework for understanding
how AD pathology may spread across cortical layers.


Pathway Enrichment Summary

L2/3 IT neurons: synaptic transmission, oxidative phosphorylation, Alzheimer disease pathway
Astrocytes: stress response, cytokine signaling, neuroinflammation, reactive gliosis
Oligodendrocytes: myelination, lipid metabolism, axon ensheathment
Microglia-PVM: immune response, phagocytosis, cytokine production
Vip interneurons: inhibitory synapse assembly, GABA signaling
Pvalb interneurons: fast synaptic transmission, potassium channel activity


REPOSITORY STRUCTURE
--------------------

seaad-mtg-scrna/
  README.md
  .gitignore
  scripts/
    load_data.R               Load h5ad file and create Seurat object
    02_normalize_pca.R        Normalize, find variable features, run PCA
    02_normalize_pca.sh       SLURM job script for normalization step
    03_umap_cluster.R         UMAP visualization and Louvain clustering
    03_umap_cluster.sh        SLURM job script for clustering step
    04_markers.R              Find marker genes per cluster
    04_markers.sh             SLURM job script for marker identification
    05_annotate.R             Apply official SEA-AD cell type labels
    05_annotate.sh            SLURM job script for annotation step
    06_de_analysis.R          Differential expression AD versus normal
    06_de_analysis.sh         SLURM job script for DE analysis
    07_pathway_enrichment.R   GO and KEGG pathway enrichment analysis
    07_pathway_enrichment.sh  SLURM job script for pathway enrichment
    08_heatmap.R              Heatmap of top DE genes
    08_heatmap.sh             SLURM job script for heatmap generation
    09_cell_proportion.R      Cell proportion analysis across groups
    09_cell_proportion.sh     SLURM job script for proportion analysis
    10_trajectory.R           Slingshot pseudotime trajectory analysis
    10_trajectory.sh          SLURM job script for trajectory analysis
  results/
    qc/
      qc_violin.png           QC metrics violin plots before filtering
      qc_scatter.png          QC scatter plots showing metric relationships
    clustering/
      elbow_plot.png          PCA elbow plot for PC selection
      all_markers.csv         All marker genes per cluster
      top5_markers.csv        Top 5 marker genes per cluster
    umap/
      umap_clusters.png       UMAP colored by cluster number 0 through 26
      umap_celltypes.png      UMAP colored by broad cell type category
      umap_cogn.png           UMAP colored by cognitive status
      umap_annotated.png      UMAP with official SEA-AD subclass labels
    differential_expression/
      plots/                  Volcano plots for each cell type
      tables/
        de_results.csv        All DE results across all cell types
        top_de_genes.csv      Top 10 DE genes per cell type
    pathway_enrichment/
      plots/                  GO and KEGG dot plots for each cell type
      tables/                 GO and KEGG CSV tables for each cell type
    heatmap/
      heatmap_top_de_genes.png
    cell_proportion/
      cell_proportion_barplot.png
      cell_proportion_dotplot.png
      cell_proportions.csv
    trajectory/
      trajectory_pseudotime.png
      trajectory_celltypes.png


HOW TO RUN
----------

Prerequisites

Access to Indiana University Quartz HPC cluster with SLURM scheduler.
R version 4.4.1 loaded via the module system.
Data downloaded to your scratch directory.

Module Setup (run before starting R or submitting any job)

module load gnu/9.3.0
module load r/4.4.1
export LD_PRELOAD=/N/soft/rhel8/gcc/14.2.0/lib64/libstdc++.so.6
export LD_LIBRARY_PATH=/N/soft/rhel8/gcc/14.2.0/lib64:$LD_LIBRARY_PATH

Download Data

The dataset is publicly available from the AWS S3 bucket:
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

Running the Pipeline

Submit each SLURM job in order and wait for each to complete before submitting the next.

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

Each SLURM script is configured to send email notifications on job start, end, and failure.

Resource Requirements Per Step

Data loading: 128 GB RAM, 8 cores, 6 hours
Normalization and PCA: 256 GB RAM, 8 cores, 12 hours
UMAP and clustering: 128 GB RAM, 8 cores, 6 hours
Marker genes: 128 GB RAM, 8 cores, 6 hours
Annotation: 64 GB RAM, 4 cores, 2 hours
Differential expression: 128 GB RAM, 8 cores, 8 hours
Pathway enrichment: 64 GB RAM, 4 cores, 4 hours
Heatmap: 64 GB RAM, 4 cores, 2 hours
Cell proportion: 64 GB RAM, 4 cores, 2 hours
Trajectory: 64 GB RAM, 4 cores, 4 hours

Monitoring Jobs

Check job status: squeue -u $USER
Check estimated start time: squeue -u $USER --start
View job log: cat /path/to/logs/jobname_JOBID.log
Cancel a job: scancel JOBID

Important Notes

The h5ad file uses AnnData format from Python and scanpy. Because zellkonverter
had GCC library compatibility issues on Quartz, we use H5SparseMatrix from HDF5Array
instead, which reads the file directly in R without any Python dependency.

The GCC mismatch between the R module compiled with gcc 9.3.0 and newer packages
requires setting LD_PRELOAD to the gcc 14.2.0 libstdc++ before starting R.
This environment variable is included in all SLURM job scripts automatically.

Scratch storage on IU Quartz has a 60-day inactivity deletion policy. Always
copy important output files such as RDS objects to home or project directories
for long-term storage once the pipeline completes.

The workshop allocation account gives better queue priority than the r00579
allocation on Quartz. Use -A workshop in SLURM job scripts for faster scheduling.


REFERENCES
----------

Leng K, et al. (2021). Molecular characterization of selectively vulnerable neurons
in Alzheimer's disease. Nature Neuroscience, 24(2), 276-287.

Allen Institute for Brain Science. (2024). Seattle Alzheimer's Disease Brain Cell
Atlas. https://sea-ad-single-cell-profiling.s3.amazonaws.com

Hao Y, et al. (2024). Dictionary learning for integrative, multimodal and scalable
single-cell analysis. Nature Biotechnology, 42, 293-304. (Seurat v5)

Street K, et al. (2018). Slingshot: cell lineage and pseudotime inference for
single-cell transcriptomics. BMC Genomics, 19, 477.

Wu T, et al. (2021). clusterProfiler 4.0: A universal enrichment tool for
interpreting omics data. The Innovation, 2(3), 100141.

Morgan M, et al. (2023). HDF5Array: HDF5 backend for DelayedArray objects.
Bioconductor.
