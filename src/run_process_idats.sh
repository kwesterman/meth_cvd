#!/bin/bash

#SBATCH -p batch,largemem
#SBATCH -c 16
#SBATCH --mem 120G
#SBATCH -t 03:00:00
#SBATCH -o "process_idats-%A.out"

module load R/3.3.2
time Rscript process_idats.R
