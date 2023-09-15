#!/bin/bash
#$ -cwd
#$ -l mem_free=30G,h_vmem=30G,h_fsize=100G
#$ -N discordance_diag
#$ -o logs_diagnostics/discord_diag.txt
#$ -e logs_diagnostics/discord_diag.txt
#$ -m e
#$ -t 5-20
echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"
## List current modules for reproducibility
module list
module load conda_R/devel
Rscript discordance_diagnostics_bayesspace_harmony.R
