#!/bin/bash

#SBATCH -p batch,largemem
#SBATCH --mem 60G
#SBATCH -c 32
#SBATCH -t 10:00:00
#SBATCH -o module_ewas-%A.out

module load R/3.3.2
time Rscript -e "Sys.setenv(RSTUDIO_PANDOC='/opt/shared/RStudio/0.98/bin/pandoc'); rmarkdown::render('module_ewas.Rmd', output_dir='../output/')"
