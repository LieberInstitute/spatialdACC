#!/bin/bash
#SBATCH --job-name=PRECAST
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=30G
#SBATCH --output=logs_precast/vis_nnsvg_batch_k.%a.txt
#SBATCH --error=logs_precast/vis_nnsvg_batch_k.%a.txt
#SBATCH --array=5-20
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

## Load the R module (absent since the JHPCE upgrade to CentOS v7)
module load conda_R

## List current modules for reproducibility
module list

## Edit with your job command
Rscript vis_precast.R

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
