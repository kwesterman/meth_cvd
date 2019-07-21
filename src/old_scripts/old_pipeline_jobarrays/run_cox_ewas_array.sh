#!/bin/bash

#SBATCH --array=1-30
#SBATCH -p batch
#SBATCH -c 4
#SBATCH --mem 20G
#SBATCH -t 01:00:00

module load R/3.3.2
time Rscript cox_ewas_array.R $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_MAX 4
