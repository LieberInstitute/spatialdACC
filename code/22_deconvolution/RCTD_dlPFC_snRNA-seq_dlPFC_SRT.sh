#!/bin/bash
#SBATCH --job-name=rctd
#SBATCH --ntasks=1                 # Run on a single CPU
#SBATCH --mem=200G
#SBATCH --cpus-per-task=30         # Number of CPU cores per task
#SBATCH --mail-type=END
#SBATCH --mail-user=kinnaryshahh@gmail.com
#SBATCH --time=1-00:00:00

echo "**** Job starts ****"
date
echo "**** SLURM info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load conda_R/4.5.x

## List current modules for reproducibility
module list

## Edit with your job command
Rscript RCTD_dlPFC_snRNA-seq_dlPFC_SRT.R

echo "**** Job ends ****"
date

