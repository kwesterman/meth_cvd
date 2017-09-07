output/meth_data_prep.pdf : meth_data_prep.Rmd
	Rscript meth_data_prep.Rmd

output/assemble_covariates.pdf : assemble_covariates.Rmd meth_data_prep.Rmd
	Rscript assemble_covariates.Rmd

output/cox_ewas.pdf : cox_ewas.Rmd assemble_covariates.Rmd meth_data_prep.Rmd
	Rscript cox_ewas.Rmd
