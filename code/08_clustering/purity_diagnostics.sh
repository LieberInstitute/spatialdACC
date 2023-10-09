#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=30G
#$SBATCH --output=logs_diagnostics/plots_diag.txt
#$SBATCH --error=logs_diagnostics/plots_diag.txt
#SBATCH --mail-type=END
#SBATCH --mail-user=kinnaryshahh@gmail.com

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"
## List current modules for reproducibility
module list
module load conda_R
Rscript purity_diagnostics.R
