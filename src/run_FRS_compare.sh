#!/bin/bash

#SBATCH -p batch
#SBATCH --mem 25G
#SBATCH -t 00:15:00

module load R/3.3.2
Rscript FRS_compare.R
