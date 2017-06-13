#!/bin/bash

#SBATCH --array=1-20
#SBATCH --time=00:15:00
#SBATCH --partition=batch
#SBATCH --mem-per-cpu=10000

module load R/3.3.2
Rscript data_readIn.R $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_MAX

# for idx in $(seq 5)
# do
# 	Rscript data_readIn.R $idx 5
# done
