Exploration of relationships between DNA methylation and cardiovascular disease (CVD) risk. Consists of two main sub-projects: 
1. An expanded epigenome-wide association study incorpoating region- and module-based analyses (described in Westerman et al. *Clin. Epigenetics* 2019; doi: https://doi.org/10.1186/s13148-019-0705-2)
2. The development and testing of a methylation-based risk score (MRS) using a penalized regression approach (described in Westerman et al. *J. Am. Heart Assoc.* 2020; doi: https://doi.org/10.1161/JAHA.119.015299)

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

**To calculate the MRS described in Westerman et al. 2020**: Use the R script and associated CpG weights files in output/mrs_calculation/.
