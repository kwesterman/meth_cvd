# Code and instructions for assessing the performance of the methylation risk
# score for cardiovascular disease in the REGICOR nested case-control for MI
#
# Datasets required:
# - "betas": N_CpGs x N_samples matrix 
#   (row names = cg IDs, column names = study-specific subject IDs)
# - "metadata": N_samples x N_covariates data frame containing technical and
#   biological covariates

library(tidyverse)
library(pROC)

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

test_mrs <- function(meta, mrs, model_spec_string, return_fit=F) {
  # Perform a logistic regression to assess the discriminative power of the MRS
  # Input: meta (metadata data frame), mrs (risk score vector),
  #   model_spec_string (string representing the model formula)
  # Output: coefficient estimates from the linear model fit
  meta$mrs <- scale(mrs)  # Standardizes MRS before regression
  glm_fit <- glm(as.formula(model_spec_string), data=meta, family="binomial")
  glm_res <- summary(glm_fit)$coef
  if (return_fit) glm_fit else glm_res
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

calc_FRS <- function(pData) {
  # Calculation of Framingham risk score based on D'Agostino 2008 
  #   doi: https://doi.org/10.1161/CIRCULATIONAHA.107.699579
  # Input: data frame with columns sex, age, chol, hdl, 
  #   ht_med (hypertension medication T/F), sbp (sys. blood pressure), 
  #   smoking (current smoker T/F), diabetes (T/F)
  # Output: vector of Framingham Risk Scores (10-year risk of CVD)
  
  FRS_data <- pData %>%
    mutate(age=pmin(pmax(age, 30), 74),  # Constrain continuous values to within specific ranges
           chol=pmin(pmax(chol, 100), 405),
           hdl=pmin(pmax(hdl, 10), 100),
           sbp=pmin(pmax(sbp, 90), 200)) %>%
    mutate(logAge=log(age),  # Take logs of continuous values
           logTC=log(chol),
           logHDL=log(hdl),
           logSBP=log(sbp),
           smoking=smk_now) %>%
    select(sex, logAge, logTC, logHDL, ht_med, logSBP, smoking, diabetes)
  
  with(FRS_data, {
    weightedSum <- ifelse(
      sex=="M",
      (3.06117 * logAge + 1.12370 * logTC - 0.93263 * logHDL + 
         1.93303 * logSBP * (1 - ht_med) + 1.99881 * logSBP * (ht_med) + 
         0.65451 * smoking + 0.57367 * diabetes),
      (2.32888 * logAge + 1.20904 * logTC - 0.70833 * logHDL + 
         2.76157 * logSBP * (1 - ht_med) + 2.82263 * logSBP * (ht_med) + 
         0.52973 * smoking + 0.69154 * diabetes))
    frs <- ifelse(
      sex=="M",
      1 - 0.88936 ^ (exp(weightedSum - 23.9802)),
      1 - 0.95012 ^ (exp(weightedSum - 26.1931)))
    frs
  })
}

# Setup and MRS calculation
stopifnot(all(metadata$subject_id == colnames(betas)))
ssl_coefs <- readRDS("ssl_mrs_coefs.rds")  # List of CpG weight vectors per study
study_weights <- readRDS("csl_weights.rds")  # Named vector of study weights
metadata$mrs <- assemble_csl(metadata, betas, ssl_coefs, study_weights)$mrs
combined_model_coefs <- readRDS("combined_mrs_coefs.rds")  # Named vector of coefficient weights
metadata$mrs_combined <- calc_mrs(combined_model_coefs, betas)
metadata$frs <- calc_FRS(metadata)  # Requires certain variables -- see calc_FRS func.

# Initial evaluation of MRS performance
# (Note: assumes case/control variable is called "case" -- alter as needed)
# Series of models to test using test_mrs as above:
model_strings <- list(
  unadjusted="case ~ mrs",
  # basic="case ~ mrs + age + sex + cell counts",  # Add cell count vars here
  # plus_risk_factors=paste0(
  #   "case ~ mrs + age + sex + cell counts + ",  # Add cell count vars here
  #   "bmi + diabetes + smk_now + ldl + hdl + sbp"),
  FRS_only="case ~ mrs + frs")
model_results <- lapply(model_strings, function(ms) test_mrs(metadata, mrs, ms))
combined_model_results <- lapply(
  model_strings, function(ms) test_mrs(metadata, mrs_combined, ms))
save("csl_model_results", "combined_model_results",
     file="main_mrs_results.RData")

# ROC curves
metadata_roc <- filter(metadata, !is.na(frs))
mrs_only_model <- test_mrs(metadata_roc, metadata_roc$mrs, 
                           "case ~ mrs", return_fit=T)
mrs_only_ROC <- roc(metadata_roc$case, predict(mrs_only_model, type="response"))
frs_model <- glm(case ~ frs, data=metadata_roc, family="binomial")
frs_ROC <- roc(metadata_roc$case, predict(frs_model, type="response"))
combined_model <- test_mrs(metadata_roc, metadata_roc$mrs, 
                           "case ~ mrs + frs", return_fit=T)
combined_ROC <- roc(metadata_roc$case, predict(combined_model, type="response"))
delong_test <- roc.test(frs_ROC, combined_ROC)
saveRDS(delong_test, "delong_test.rds")

png("roc_curves.png")
plot(mrs_only_ROC, col="blue")
plot(frs_ROC, col="red", add=T)
plot(combined_ROC, col="green", add=T)
legend("bottomright", 
       legend=paste0(c("MRS: auc=", "FRS: auc=", "MRS+FRS: auc="), 
                     round(as.numeric(c(mrs_only_ROC$auc, frs_ROC$auc, 
                                        combined_ROC$auc)), 2)), 
       col=c("blue", "red", "green"), lty=1, lwd=3, cex=1.3)
dev.off()

# Age, sex, and Framingham Risk Score interactions
test_mrs_with_SE <- function(nmd, mrs, subset, covar_string) {
  nmd$mrs <- mrs
  subset_data <- nmd[subset, ]
  test_res <- test_mrs(subset_data, subset_data$mrs, 
                       as.formula(paste("case ~ mrs +", covar_string)),
                       return_fit=T)
  model_coefs <- summary(test_res)$coef["mrs", c("Estimate", "Std. Error")]
  tibble(OR_per_SD=exp(model_coefs["Estimate"]),
         SE_lower=exp(model_coefs["Estimate"] - model_coefs["Std. Error"]),
         SE_higher=exp(model_coefs["Estimate"] + model_coefs["Std. Error"]))
}

sex_strat <- map_dfr(unique(metadata$sex), function(s) {
  tryCatch(test_mrs_with_SE(metadata, mrs, metadata$sex == s, 
                            "age"),
           error=function(e) data.frame(OR_per_SD=NA, SE_lower=NA, SE_higher=NA))
}, .id="Sex") %>%
  mutate(Sex=unique(metadata$sex))  # Ensure sex label order is correct

metadata$age_group <- cut(metadata$age,
                          breaks=c(0, 55, 65, 75, 100),
                          labels=c("0-55", "56-65", "66-75", "76+"))
age_group_strat <- map_dfr(levels(metadata$age_group), function(g) {
  tryCatch(test_mrs_with_SE(metadata, mrs, metadata$age_group == g,
                            "age + sex"),
           error=function(e) data.frame(OR_per_SD=NA, SE_lower=NA, SE_higher=NA))
}, .id="Age group") %>%
  mutate(`Age group`=levels(metadata$age_group))

metadata$frs_group <- cut(metadata$frs,
                          breaks=c(0, 0.1, 0.2, 0.4, 1),
                          labels=c("0-0.1", "0.1-0.2", "0.2-0.4", "0.4+"))
frs_group_strat <- map_dfr(levels(metadata$frs_group), function(g) {
  tryCatch(test_mrs_with_SE(metadata, mrs, metadata$frs_group == g,
                            "age + sex"),
           error=function(e) data.frame(OR_per_SD=NA, SE_lower=NA, SE_higher=NA))
}, .id="FRS risk group") %>%
  mutate(`FRS risk group`=levels(metadata$frs_group))

save("sex_strat", "age_group_strat", "frs_group_strat",
     file="stratified_results.RData")

