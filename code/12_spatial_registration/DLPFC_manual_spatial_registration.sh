#!/bin/bash
#SBATCH --job-name=spatial_reg
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --output=logs/DLPFC_spatial_reg_k.%a.txt
#SBATCH --error=logs/DLPFC_spatial_reg_k.%a.txt
#SBATCH --array=7-10
#SBATCH --mail-type=END
#SBATCH --mail-user=kinnaryshahh@gmail.com

echo "**** Job starts ****"
date
echo "**** SLURM info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load conda_R

## List current modules for reproducibility
module list

## Edit with your job command
Rscript DLPFC_manual_spatial_registration.R

echo "**** Job ends ****"
date
