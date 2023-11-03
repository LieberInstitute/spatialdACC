#!/bin/bash
#SBATCH --job-name=bayesspace_mnn
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=30G
#SBATCH --array=5-20
#SBATCH --mail-type=END
#SBATCH --mail-user=kinnaryshahh@gmail.com

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## Load the R module (absent since the JHPCE upgrade to CentOS v7)
module load conda_R

## List current modules for reproducibility
module list

## Edit with your job command
Rscript bayesspace_mnn.R

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
