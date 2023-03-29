#!/bin/bash
#$ -cwd
#$ -l mem_free=100G,h_vmem=100G,h_stack=256M,h_fsize=100G
#$ -o /dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/code/VistoSeg/code/logs/$TASK_ID_splitSlide.txt 
#$ -e /dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/code/VistoSeg/code/logs/$TASK_ID_splitSlide.txt
#$ -m e
#$ -M ryan.miller@libd.org
#$ -t 5
#$ -tc 4

echo "**** Job starts ****"
date


echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"s
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## load MATLAB
module load matlab/R2019a

## Load toolbox for VistoSeg
toolbox='/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/code/VistoSeg/code/'

## Read inputs from splitSlide_list.txt file
fname=$(awk 'BEGIN {FS="\t"} {print $1}' splitSlide_list.txt | awk "NR==${SGE_TASK_ID}")
A1=$(awk 'BEGIN {FS="\t"} {print $2}' splitSlide_list.txt | awk "NR==${SGE_TASK_ID}")
B1=$(awk 'BEGIN {FS="\t"} {print $3}' splitSlide_list.txt | awk "NR==${SGE_TASK_ID}")
C1=$(awk 'BEGIN {FS="\t"} {print $4}' splitSlide_list.txt | awk "NR==${SGE_TASK_ID}")
D1=$(awk 'BEGIN {FS="\t"} {print $5}' splitSlide_list.txt | awk "NR==${SGE_TASK_ID}")

## Run refineVNS function
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), splitSlide('$fname',$A1,$B1,$C1,$D1)"

echo "**** Job ends ****"
date



