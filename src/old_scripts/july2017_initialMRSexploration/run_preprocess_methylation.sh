#!/bin/bash

#SBATCH -p batch,largemem
#SBATCH -c 8
#SBATCH --mem 120G
#SBATCH -t 4:00:00
#SBATCH -o preprocess_methylation_%A.out

module load R/3.3.2
time Rscript preprocess_methylation.R
