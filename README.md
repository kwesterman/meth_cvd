Development of a methylation-based risk score, and exploration of its performance in relation to other "dimensions" of cardiovascular risk based on demographics, clinical risk factors, and predicted genetic risk.

metadata: Reads in demographic and phenotypic data from a wide variety of sources, followed by data cleaning and integration into straightforward metadata datasets for incorporation with methylation and genotypes.

process_idats: Reads in raw .idat files, performs array-wide signal detection tests, and produces an intensity-based sample QC plot.

methylation_preprocessing: Performs sample QC, Noob preprocessing (background and dye bias correction), and BMIQ within-sample normalization. Main output is the "final" matrix of M-values to be used in downstream analyses.

substructure: Calculation of various quantities relevant to examining and adjusting for substructure in the methylation dataset. Includes cell count estimation, PCA calculation, CPACOR (PCA on control probes), and ComBat.

module_ewas: RMarkdown script for generating the results/manuscript.
