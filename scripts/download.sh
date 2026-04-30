#!/bin/bash
#SBATCH --job-name=download_seaad
#SBATCH -A r00579
#SBATCH --partition=general
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=08:00:00
#SBATCH --output=/N/scratch/dikonda/seurat_project/logs/download_%j.log

cd /N/scratch/dikonda/seurat_project/data
wget https://sea-ad-single-cell-profiling.s3.amazonaws.com/MTG/RNAseq/SEAAD_MTG_RNAseq_all-nuclei.2024-02-13.h5ad
