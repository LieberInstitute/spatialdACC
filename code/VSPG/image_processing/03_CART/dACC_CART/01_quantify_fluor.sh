#!/bin/bash
#SBATCH --mem=20G
#SBATCH --job-name=quantify_fluo
#SBATCH --array=1-4%4
#SBATCH -o /dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/code/VSPG/image_processing/03_CART/dACC_CART/logs/1_quant_fluo%a.log
#SBATCH -e /dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/code/VSPG/image_processing/03_CART/dACC_CART/logs/1_quant_fluo%a.log


echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

module load cellpose/2.0
python 01_quantify_fluor.py

echo "**** Job ends ****"
date