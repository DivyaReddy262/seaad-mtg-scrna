#!/bin/bash
#SBATCH --job-name=cell_prop
#SBATCH -A workshop
#SBATCH --partition=general
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --output=/N/scratch/dikonda/seurat_project/logs/cell_prop_%j.log
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=dikonda@iu.edu

module load gnu/9.3.0
module load r/4.4.1
export LD_PRELOAD=/N/soft/rhel8/gcc/14.2.0/lib64/libstdc++.so.6
export LD_LIBRARY_PATH=/N/soft/rhel8/gcc/14.2.0/lib64:$LD_LIBRARY_PATH

Rscript /N/scratch/dikonda/seurat_project/scripts/09_cell_proportion.R
