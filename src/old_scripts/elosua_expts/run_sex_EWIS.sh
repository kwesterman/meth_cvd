#!/bin/bash

#SBATCH -p batch,largemem
#SBATCH --mem 40G
#SBATCH -c 32
#SBATCH -t 1:00:00
#SBATCH -o sex_EWIS-%A.out

module load R/3.4.3
time Rscript -e "Sys.setenv(RSTUDIO_PANDOC='/opt/shared/RStudio/0.98/bin/pandoc'); rmarkdown::render('sex_EWIS.Rmd', output_dir='../output/')"
