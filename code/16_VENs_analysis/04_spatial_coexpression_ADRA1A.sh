#!/bin/bash
#SBATCH --job-name=coexpr
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=25G
#SBATCH --mail-type=END
#SBATCH --output=logs/coexpr_ADRA1A.txt
#SBATCH --error=logs/coexpr_ADRA1A.txt
#SBATCH --mail-user=kinnaryshahh@gmail.com
#SBATCH --time=3-00:00:00

module load conda_R/4.3.x
Rscript 04_spatial_coexpression_ADRA1A.R
