#!/bin/bash
#$ -cwd
#$ -l mem_free=30,h_vmem=30,h_fsize=100G
#$ -pe local 12
#$ -N nnSVG
#$ -o logs_nnSVG/nnSVG_1000.txt
#$ -e logs_nnSVG/nnSVG_1000.txt
#$ -m e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

module load conda_R/devel

## List current modules for reproducibility
module list

## Edit with your job command
Rscript nnSVG.R

echo "**** Job ends ****"
date
