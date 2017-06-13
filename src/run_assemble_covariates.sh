#!/bin/bash

#SBATCH -p batch,largemem
#SBATCH --mem 100G
#SBATCH -t 2:00:00
#SBATCH -o assemble_covariates_%A.out

module load R/3.3.2
time Rscript assemble_covariates.R
