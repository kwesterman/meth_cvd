#!/bin/bash

#SBATCH -p batch,largemem
#SBATCH --mem 50G
#SBATCH -n 32
#SBATCH -t 2:00:00
#SBATCH -o past_ewas-%A.out

module load R/3.3.2
time Rscript -e "Sys.setenv(RSTUDIO_PANDOC='/opt/shared/RStudio/0.98/bin/pandoc'); rmarkdown::render('past_ewas.Rmd, output_dir='../output/')"
