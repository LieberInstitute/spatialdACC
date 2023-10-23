#!/bin/bash
#SBATCH --job-name=nnSVG_PRECAST
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --output=logs/nnsvg_PRECAST_DE_k.%a.txt
#SBATCH --error=logs/nnsvg_PRECAST_DE_k.%a.txt
#SBATCH --array=5-12
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
Rscript nnSVG_precast_DE.R

echo "**** Job ends ****"
date
