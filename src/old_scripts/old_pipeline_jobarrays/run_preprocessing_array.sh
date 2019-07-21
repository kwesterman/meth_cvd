#!/bin/bash

#SBATCH --array=1-20
#SBATCH --partition=batch
#SBATCH -c 4
#SBATCH --mem=40G
#SBATCH --time=01:00:00

module load R/3.3.2
Rscript preprocess_methylation.R $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_MAX 4

