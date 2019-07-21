#!/bin/bash

#SBATCH -p batch,largemem
#SBATCH --mem 90G
#SBATCH -c 1
#SBATCH -t 4:00:00
#SBATCH -o mrs_experiments-%A.out

module load R/3.3.2
time Rscript -e "Sys.setenv(RSTUDIO_PANDOC='/opt/shared/RStudio/0.98/bin/pandoc'); rmarkdown::render('mrs_experiments.Rmd', output_dir='../output/')"
