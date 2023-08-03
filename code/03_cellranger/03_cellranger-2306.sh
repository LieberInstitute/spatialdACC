#!/bin/bash
#$ -cwd
#$ -pe local 8
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -N spatialdACC_cellranger
#$ -o logs/cellranger-2306.1.txt
#$ -e logs/cellranger-2306.1.txt
#$ -t 1
#$ -tc 1

# Files to process
# 5c_dACC_MRV 

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
SAMPLE=$(awk 'BEGIN {FS="\t"} {print $1}' samples.txt | awk "NR==${SGE_TASK_ID}")
SAMPLEID=$(awk 'BEGIN {FS="\t"} {print $2}' samples.txt | awk "NR==${SGE_TASK_ID}")
## SAMPLE=$(awk "NR==${SGE_TASK_ID}" 03_cellranger.txt)
echo "Processing sample ${SAMPLE}"
date

## Run CellRanger
cellranger count --id=${SAMPLE} \
    --transcriptome=/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A \
    --fastqs=/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/raw-data/FASTQ/2023-06-09_SPag051623/${SAMPLE} \
    --sample=${SAMPLEID} \
    --jobmode=local \
    --localcores=8 \
    --localmem=64

## Move output
echo "Moving data to new location"
date
mkdir -p /dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/processed-data/03_cellranger/${SAMPLE}
mv ${SAMPLE} /dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/processed-data/03_cellranger/

echo "**** Job ends ****"
date
