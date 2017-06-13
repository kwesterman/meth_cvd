#!/bin/bash

#SBATCH -p batch
#SBATCH -c 16
#SBATCH --mem 50G
#SBATCH -t 02:00:00
#SBATCH -o "process_idats_%A.out"

module load R/3.3.2
time Rscript process_idats.R 16
