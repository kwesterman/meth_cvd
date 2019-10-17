Exploration of relationships between DNA methylation and cardiovascular disease (CVD) risk. Consists of two main sub-projects: an expanded epigenome-wide association study incorpoating region- and module-based analyses, and the development of a methylation-based risk score (MRS).

src/ contains analysis scripts.
output/ contains output from other programs and individuals.
doc/ contains manuscript files and outputs.

* 0_metadata: Reads in demographic and phenotypic data from a wide variety of sources, followed by data cleaning and integration into straightforward metadata datasets for incorporation with methylation and genotypes.
* 1_process_idats: Reads in raw .idat files, performs array-wide signal detection tests, and produces an intensity-based sample QC plot.
* 2_methylation_preprocessing: Performs sample QC, Noob preprocessing (background and dye bias correction), and BMIQ within-sample normalization. Main output is the "final" matrix of M-values to be used in downstream analyses.
* 3_substructure: Calculation of various quantities relevant to examining and adjusting for substructure in the methylation dataset. Includes cell count estimation, PCA calculation, CPACOR (PCA on control probes), and ComBat.
* clean_genotypes_*.sh and calc_grs.sh: Genotype preprocessing and CAD genetic risk score calculation for integration with MRS.
* module_ewas and module_ewas_supp: RMarkdown scripts for generating the results and supplementary information for the module-based incident CVD EWAS manuscript.
* mrs_models (in doc/mrs_models/): RMarkdown script for generating the main results for the MRS manuscript.
