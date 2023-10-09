#!/bin/bash
#!/bin/bash
#SBATCH --job-name=nnSVG_PRECAST
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --output=logs_diagnostics/discord_diag_nnSVG_precast_k.%a.txt
#SBATCH --error=logs_diagnostics/discord_diag_nnSVG_precast_k.%a.txt
#SBATCH --array=5-20
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
Rscript discordance_diagnostics_nnSVG_precast.R
