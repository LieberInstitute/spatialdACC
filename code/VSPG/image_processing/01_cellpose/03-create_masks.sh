#!/bin/bash
#SBATCH --reservation=caracol2
#SBATCH --partition=caracol
#SBATCH --gpus=1
#SBATCH --mem=200G
#SBATCH --job-name=create_mask
#SBATCH -o /dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/code/VSPG/image_processing/logs/1-create_masks_%a.log
#SBATCH -e /dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/code/VSPG/image_processing/logs/1-create_masks_%a.log
#SBATCH --array=2-4%2

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

###############################################################################
#   Dynamically select a GPU based on availability
###############################################################################

#USAGE_CUTOFF=10
#NUM_GPUS=1
#
#avail_gpus=$(nvidia-smi --query-gpu=utilization.gpu --format=csv,noheader | cut -d " " -f 1 | awk -v usage="$USAGE_CUTOFF" '$1 < usage {print NR - 1}')
#
##  Simply exit with an error if there are no GPUs left
#if [[ -z $avail_gpus ]]; then
#    echo "No GPUs are available."
#    exit 1
#fi
#
#export CUDA_VISIBLE_DEVICES=$(echo "$avail_gpus" | head -n $NUM_GPUS | paste -sd ",")
#
#echo "Chose GPU(s): $CUDA_VISIBLE_DEVICES"
#
###############################################################################
#   Submit the python script
###############################################################################

module load cellpose/2.2

#tif_dir="/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/processed-data/Images/VistoSeg/Capture_areas/if-images/"
#tif_files=($tif_dir*.tif) 
#sample_path=${tif_files[$SLURM_ARRAY_TASK_ID-1]}
#sample_id=$(basename "$sample_path" .tif)
#echo $sample_id
echo "starting python"
python /dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/code/VSPG/image_processing/01_cellpose/03-create_masks.py

echo "**** Job ends ****"
date

