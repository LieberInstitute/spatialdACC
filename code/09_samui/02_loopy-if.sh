#!/bin/bash
#SBATCH --mem=10G
#SBATCH --job-name=02samui
#SBATCH --output=logs/%x.txt
#SBATCH --error=logs/%x.txt
# #SBATCH --array=1-4

#echo -e "\n"
echo "date: $(date)"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURM_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load samui/1.0.0-next.45
python 02_loopy-if.py

echo "**** Job ends ****"
date

#echo -e "\n"
