#!/bin/bash
#SBATCH --job-name=nmf
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=15G
#SBATCH --output=logs/run_NMF.txt
#SBATCH --error=logs/run_NMF.txt
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
Rscript 02_run_NMF.R

echo "**** Job ends ****"
date
