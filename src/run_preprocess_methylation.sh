#!/bin/bash

#SBATCH -p batch,largemem
#SBATCH -c 8
#SBATCH --mem 180G
#SBATCH -t 6:00:00
#SBATCH -o preprocess_methylation-%A.out

module load R/3.3.2
time Rscript preprocess_methylation.R
