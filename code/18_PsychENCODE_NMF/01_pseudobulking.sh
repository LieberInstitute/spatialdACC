#!/bin/bash
#SBATCH --job-name=pseudo
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=750G
#SBATCH --output=logs/pseudo.txt
#SBATCH --error=logs/pseudo.txt
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

module load conda_R/4.4.x

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 01_pseudobulking.R

echo "**** Job ends ****"
date
