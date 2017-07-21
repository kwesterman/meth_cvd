Pertains to Aim 1 of my thesis, consisting of the construction and analysis of a methylation risk score predicting cardiovascular disease incidence. Because some steps require visual/manual intervention, the basic methylation preprocessing pipeline is run by hand -- a basic description is below:

process_idats: Reads in raw .idat files, performs array-wide signal detection tests, and produces an intensity-based sample QC plot.

methylation_preprocessing: Performs sample QC, Noob preprocessing (background and dye bias correction), and BMIQ within-sample normalization. Main output is the "final" matrix of M-values to be used in downstream analyses.

assemble_covariates: Assemble the series of covariates to be used in downstream regression analysis. Some come from cohort phenotypes data files, and some are calculated from the methylation data (e.g. PCs). Main output is the matrix of all non-methylation covariates.
