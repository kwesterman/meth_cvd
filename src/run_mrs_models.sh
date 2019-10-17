#!/bin/bash

#SBATCH -p batch,largemem
#SBATCH --mem 220G
#SBATCH -c 5
#SBATCH -t 36:00:00
#SBATCH -o mrs_models-%A.out

module load R/3.4.3
time Rscript -e "Sys.setenv(RSTUDIO_PANDOC='/opt/shared/RStudio/0.98/bin/pandoc'); rmarkdown::render('mrs_models.Rmd', output_dir='../output/')"
