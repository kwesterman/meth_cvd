#!/bin/bash

#SBATCH -p batch,largemem
#SBATCH --mem 100G
#SBATCH -t 3:00:00
#SBATCH -o substructure_ewas-%A.out

module load R/3.3.2
time Rscript substructure_ewas.R whi fhs
