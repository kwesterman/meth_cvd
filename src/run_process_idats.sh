#!/bin/bash

#SBATCH -p batch,largemem
#SBATCH -n 16
#SBATCH --mem 120G
#SBATCH -t 05:00:00
#SBATCH -o "process_idats-%A.out"

module load R/3.3.2
time Rscript process_idats.R
