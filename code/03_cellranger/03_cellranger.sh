#!/bin/bash
#$ -cwd
#$ -l mem_free=8G,h_vmem=8G,h_fsize=100G
#$ -pe local 4
#$ -N spatialdACC_cellranger
#$ -o logs/cellranger.$TASK_ID.txt
#$ -e logs/cellranger.$TASK_ID.txt
#$ -m e
#$ -t 1-4
#$ -tc 4

# Files to process
# 1c_dACC_MRV 
# 2c_dACC_MRV 
# 3c_dACC_MRV 
# 4c_dACC_MRV 
 

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## load CellRanger
module load cellranger/7.0.0

## List current modules for reproducibility
module list

## Locate file
SAMPLE=$(awk 'BEGIN {FS="\t"} {print $1}' 03_cellranger.txt | awk "NR==${SGE_TASK_ID}")
SAMPLEID=$(awk 'BEGIN {FS="\t"} {print $2}' 03_cellranger.txt | awk "NR==${SGE_TASK_ID}")
## SAMPLE=$(awk "NR==${SGE_TASK_ID}" 03_cellranger.txt)
echo "Processing sample ${SAMPLE}"
date

## Run CellRanger
cellranger count --id=${SAMPLE} \
    --transcriptome=/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A \
    --fastqs=/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/raw-data/FASTQ/2023-05-04_SPag041423/${SAMPLE} \
    --sample=${SAMPLEID} \
    --jobmode=local \
    --localcores=8 \
    --localmem=64

## Move output
echo "Moving data to new location"
date
mkdir -p /dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/processed-data/03_cellranger/
mv ${SAMPLE} /dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/processed-data/03_cellranger/

echo "**** Job ends ****"
date