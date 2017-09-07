#!/bin/bash

#SBATCH -p batch,largemem
#SBATCH -n 100
#SBATCH --mem 40G
#SBATCH -t 8:00:00
#SBATCH -o cox_ewas-%A.out

module load R/3.3.2
time Rscript cox_ewas.R
