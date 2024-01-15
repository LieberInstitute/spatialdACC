#!/bin/bash
#SBATCH --job-name=doublet
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --output=logs/doublet.txt
#SBATCH --error=logs/doublet.txt
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
Rscript 03_doublet_detection.R

echo "**** Job ends ****"
date
