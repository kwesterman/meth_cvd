#!/bin/bash

#SBATCH -p batch,largemem
#SBATCH -c 32
#SBATCH --mem 40G
#SBATCH -t 4:00:00
#SBATCH -o cox_ewas_%A.out

module load R/3.3.2
time Rscript cox_ewas.R
