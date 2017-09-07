for (pkg in c("knitr","tidyverse","rmarkdown","gtools")) {
  if(!require(pkg, character.only=T)) {
    install.packages(pkg)
  }
}
for (pkg in c("minfi","wateRmelon","RefFreeEWAS","FlowSorted.Blood.450k")) {
  if(!require(pkg, character.only=T)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite(pkg)
  }
}

render("meth_data_prep.Rmd", output_dir="../output/")
render("assemble_covariates.Rmd", output_dir="../output/")
render("cox_ewas.Rmd", output_dir="../output/")