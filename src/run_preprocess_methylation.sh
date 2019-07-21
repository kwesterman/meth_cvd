#!/bin/bash

#SBATCH -p batch,largemem
#SBATCH -c 8
#SBATCH --mem 160G
#SBATCH -t 6:00:00
#SBATCH -o preprocess_methylation-%A.out

module load R/3.3.2
time Rscript 2_preprocess_methylation.R whi fhs
