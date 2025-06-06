#!/bin/bash
#SBATCH --mem=20G
#SBATCH --job-name=cart
#SBATCH --array=1-4%4
#SBATCH -o /dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/code/VSPG/image_processing/03_CART/dACC_CART/logs/1_cart%a.log
#SBATCH -e /dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/code/VSPG/image_processing/03_CART/dACC_CART/logs/1_cart%a.log


echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

module load cellpose/2.0
python 04b_cart_test.py

echo "**** Job ends ****"
date