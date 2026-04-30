#!/bin/bash
#SBATCH --job-name=de_analysis
#SBATCH -A workshop
#SBATCH --partition=general
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --time=08:00:00
#SBATCH --output=/N/scratch/dikonda/seurat_project/logs/de_analysis_%j.log
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=dikonda@iu.edu

module load gnu/9.3.0
module load r/4.4.1
export LD_PRELOAD=/N/soft/rhel8/gcc/14.2.0/lib64/libstdc++.so.6
export LD_LIBRARY_PATH=/N/soft/rhel8/gcc/14.2.0/lib64:$LD_LIBRARY_PATH

Rscript /N/scratch/dikonda/seurat_project/scripts/06_de_analysis.R
