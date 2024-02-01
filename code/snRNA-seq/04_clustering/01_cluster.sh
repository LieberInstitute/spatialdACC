#!/bin/bash
#SBATCH --job-name=cluster
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100G
#SBATCH --output=logs/cluster.txt
#SBATCH --error=logs/cluster.txt
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
Rscript 01_cluster.R

echo "**** Job ends ****"
date
