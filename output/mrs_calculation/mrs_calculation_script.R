# Code and instructions for calculating the methylation risk scores for
# cardiovascular disease described in Westerman et al. JAHA 2020
#
# Datasets required (pre-load into R environment):
# - "betas": N_CpGs x N_samples matrix 
#   (row names = cg IDs, column names = study-specific subject IDs)
# - "metadata": N_samples x N_covariates data frame containing technical and
#   biological covariates


library(tidyverse)


calc_mrs <- function(coefs, meth_mat) {
  # Calculate the methlylation risk score for each subject
  # Input: coefs (named vector of CpG weights), meth_mat (beta_value matrix)
  # Output: vector of beta-value weighted sums
  available_cpgs <- names(coefs)[names(coefs) %in% rownames(meth_mat)]
  print(paste(length(available_cpgs), "of", length(coefs), 
              "MRS CpGs available."))
  mrs <- as.vector(
    t(meth_mat[available_cpgs, , drop=F]) %*% coefs[available_cpgs])
  mrs
}

csl_predict <- function(stacked_df, studies, weights) {
  # Calculate a final CSL prediction given study-specific predictions
  # Input: stacked_df (data frame of study-specific predictions),
  #   studies (list of studies), weights (weights for linear combination)
  # Output: vector of final CSL predictions
  for (s in studies) stacked_df[[s]] <- scale(stacked_df[[s]])
  as.vector(as.matrix(stacked_df[, studies]) %*% weights)
}

assemble_csl <- function(meta, betas, cpg_weights_list, csl_weights) {
  # Use a list of CpG weight vectors (one element per study)
  # to calculate SSL predictions, find study-specific weights,
  # and return a final CSL prediction
  # Calculate SSL predictions and generate final CSL prediction
  # Input: meta (metadata data frame), betas (N_CpGs x N_samples matrix),
  #   cpg_weights_list (named list of study-specific CpG coefficients),
  #   csl_weights (vector of weights for combination of SSL predictions)
  # Output: list containing the final prediction and the CSL weights used
  ssl_preds <- meta
  for (s in names(cpg_weights_list)) {
    ssl_preds[[s]] <- calc_mrs(cpg_weights_list[[s]], betas)
  }
  mrs <- csl_predict(ssl_preds, names(csl_weights), csl_weights)
  list(mrs=mrs, csl_weights=csl_weights)
}

# Setup and MRS calculation
stopifnot(all(metadata$subject_id == colnames(betas)))
ssl_coefs <- readRDS("ssl_mrs_coefs.rds")  # List of CpG weight vectors per study
study_weights <- readRDS("csl_weights.rds")  # Named vector of study weights
metadata$mrs <- assemble_csl(metadata, betas, ssl_coefs, study_weights)$mrs
combined_model_coefs <- readRDS("combined_mrs_coefs.rds")  # Named vector of coefficient weights
metadata$mrs_combined <- calc_mrs(combined_model_coefs, betas)
combat_model_coefs <- readRDS("combat_mrs_coefs.rds")
metadata$mrs_combat <- calc_mrs(combat_model_coefs, betas)
