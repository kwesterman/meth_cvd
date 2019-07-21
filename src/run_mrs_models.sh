#!/bin/bash

#SBATCH -p batch,largemem
#SBATCH --mem 180G
#SBATCH -c 8
#SBATCH -t 10:00:00
#SBATCH -o mrs_models-%A.out

module load R/3.4.3
time Rscript -e "Sys.setenv(RSTUDIO_PANDOC='/opt/shared/RStudio/0.98/bin/pandoc'); rmarkdown::render('mrs_models.Rmd', output_dir='../output/')"
