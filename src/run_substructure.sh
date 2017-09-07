#!/bin/bash

#SBATCH -p batch,largemem
#SBATCH --mem 120G
#SBATCH -t 6:00:00
#SBATCH -o substructure-%A.out

module load R/3.3.2
time Rscript substructure.R
