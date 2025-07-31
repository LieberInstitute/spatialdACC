#!/bin/bash
#SBATCH -p gpu
#SBATCH --mem=20G
#SBATCH --job-name=create_samui_postCART
#SBATCH -o logs/create_samui_dlpfc.txt 
#SBATCH -e logs/create_samui_dlpfc.txt 
#SBATCH --array=1-4%2

echo "**** Job starts ****"
date


echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"s
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURM_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load samui/1.0.0-next.24
python 03_build_samui_postCART.py 

echo "**** Job ends ****"
date