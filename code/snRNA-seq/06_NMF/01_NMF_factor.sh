#!/bin/bash
#SBATCH --job-name=nmf_factor
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=10G
#SBATCH --output=logs/project.txt
#SBATCH --error=logs/project.txt
#SBATCH --mail-type=END
#SBATCH --mail-user=kinnaryshahh@gmail.com
#SBATCH --time=2-00:00:00

echo "**** Job starts ****"
date
echo "**** SLURM info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load conda_R/4.3.x

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 01_NMF_factor.R

echo "**** Job ends ****"
date
