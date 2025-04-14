#!/bin/bash
#SBATCH --job-name=cv
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=40G
#SBATCH --output=logs/cv.txt
#SBATCH --error=logs/cv.txt
#SBATCH --mail-type=END
#SBATCH --mail-user=kinnaryshahh@gmail.com
#SBATCH --time=7-00:00:00

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
Rscript 01_cross_validation.R

echo "**** Job ends ****"
date
