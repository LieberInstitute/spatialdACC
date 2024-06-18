#!/bin/bash
#SBATCH --job-name=coexpr
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=25G
#SBATCH --mail-type=END
#SBATCH --output=logs/coexpr_DLPFC_30.txt
#SBATCH --error=logs/coexpr_DLPFC_30.txt
#SBATCH --mail-user=kinnaryshahh@gmail.com
#SBATCH --time=3-00:00:00

module load conda_R
Rscript 06_spatial_coexpression_DLPFC_30.R
