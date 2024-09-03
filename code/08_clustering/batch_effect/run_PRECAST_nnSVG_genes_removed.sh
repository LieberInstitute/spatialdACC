#!/bin/bash
#SBATCH --job-name=nnSVG_PRECAST
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=30G
#SBATCH --output=logs_precast/nnsvg_PRECAST_batch_k.%a.txt
#SBATCH --error=logs_precast/nnsvg_PRECAST_batch_k.%a.txt
#SBATCH --array=7-9
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
Rscript run_PRECAST_nnSVG_genes_removed.R

echo "**** Job ends ****"
date
