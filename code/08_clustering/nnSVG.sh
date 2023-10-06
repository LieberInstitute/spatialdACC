#!/bin/bash
#SBATCH --job-name=nnSVG
#SBATCH --output=logs_nnSVG/nnSVG.txt
#SBATCH --error=logs_nnSVG/nnSVG.txt
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=20G
#SBATCH --mail-type=END
#SBATCH --mail-user=kinnaryshahh@gmail.com # Please replace with the appropriate email address
echo "**** Job starts ****"
date
echo "**** SLURM info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load conda_R/devel

## List current modules for reproducibility
module list

## Edit with your job command
Rscript nnSVG.R

echo "**** Job ends ****"
date
