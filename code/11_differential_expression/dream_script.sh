#!/bin/bash
#SBATCH --job-name=dream
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=10G
#SBATCH --output=logs/dream.txt
#SBATCH --error=logs/dream.txt
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

module load conda_R/4.3.x

## List current modules for reproducibility
module list

## Edit with your job command
Rscript dream_script_consistent_logFC.R

echo "**** Job ends ****"
date
