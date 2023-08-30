#!/bin/bash
#SBATCH --mem=50G
#SBATCH --job-name=r2python
#SBATCH --array=1



echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURM_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load conda_R/4.3
Rscript 01_r2py.R

mkdir -p ./logs/

cat ./slurm-${SLURM_JOBID}.out >> ./logs/01_r2py.out
#rm ./slurm-${SLURM_JOBID}.out

echo "**** Job ends ****"
date