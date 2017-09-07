#!/bin/bash

#SBATCH -p batch
#SBATCH --mem 5G

module load R/3.3.2
Rscript assemble_covariates.R
