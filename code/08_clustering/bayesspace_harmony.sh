#!/bin/bash
#$ -cwd
#$ -l mem_free=30G,h_vmem=30G,h_fsize=30G
#$ -N bayesSpace_captureArea_k_many
#$ -o logs_bayesspace_harmony/bayesSpace_captureArea_k.$TASK_ID.txt
#$ -e logs_bayesspace_harmony/bayesSpace_captureArea_k.$TASK_ID.txt
#$ -m e
#$ -t 10-20
#$ -tc 2

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## Load the R module (absent since the JHPCE upgrade to CentOS v7)
module load conda_R/devel

## List current modules for reproducibility
module list

## Edit with your job command
Rscript bayesspace_harmony.R

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
