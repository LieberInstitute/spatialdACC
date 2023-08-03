#!/bin/bash
#SBATCH --mem=30G
#SBATCH --job-name=spatialdACC_build_spe
#SBATCH -o logs/raw_spe_wIFo.txt
#SBATCH -e logs/raw_spe_wIFe.txt


echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## Load the R module (absent since the JHPCE upgrade to CentOS v7)
module load conda_R/4.3

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 01_raw_spe.R

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/