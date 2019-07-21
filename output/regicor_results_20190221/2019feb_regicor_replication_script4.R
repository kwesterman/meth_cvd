# Code and instructions for assessing the performance of the methylation risk
# score for cardiovascular disease in the REGICOR nested case-control for MI
#
# Datasets required:
# - "betas": N_CpGs x N_samples matrix 
#   (row names = cg IDs, column names = study-specific subject IDs)
# - "metadata": N_samples x N_covariates data frame containing technical and
#   biological covariates

rm(list=ls())

library(tidyverse)
library(pROC)

setwd("/projects/regicor/METIAM/EPIC/kenny_replication/2019feb")

#AFS##Load betas
load("../betas.RData")

#AFS##Load metadata
load("metadata.RData")


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
#AFS# [1] "299 of 320 MRS CpGs available."
#AFS# [1] "108 of 116 MRS CpGs available."
#AFS# [1] "476 of 515 MRS CpGs available."
#AFS# [1] "340 of 369 MRS CpGs available."
combined_model_coefs <- readRDS("combined_mrs_coefs.rds")  # Named vector of coefficient weights
metadata$mrs_combined <- calc_mrs(combined_model_coefs, betas)
# [1] "93034 of 100000 MRS CpGs available."
# metadata$frs <- calc_FRS(metadata)  # Requires certain variables -- see calc_FRS func.

# Initial evaluation of MRS performance
# (Note: assumes case/control variable is called "case" -- alter as needed)
# Series of models to test using test_mrs as above:
model_strings <- list(
  unadjusted="case ~ mrs + sva1 + sva2",
  basic="case ~ mrs + age + sex + CD4T + CD8T + NK + Bcell + Mono + sva1 + sva2",
  # plus_risk_factors=paste0(
  #   "case ~ mrs + age + sex + CD4T + CD8T + NK + Bcell + Mono + ",
  #   "bmi + diabetes + smk_now + ldl + hdl + sbp")
  plus_risk_factors=paste0(
    "case ~ mrs + age + sex + CD4T + CD8T + NK + Bcell + Mono + ",
    "bmi + diabetes + smk_now + h_lipid + htn + sva1 + sva2")
  # FRS_only="case ~ mrs + frs"
)
csl_model_results <- lapply(model_strings, function(ms) {
  test_mrs(metadata, metadata$mrs, ms)
})
combined_model_results <- lapply(
  model_strings, function(ms) test_mrs(metadata, metadata$mrs_combined, ms))
save("csl_model_results", "combined_model_results",
     file="main_mrs_results.RData")
csl_model_results
# $unadjusted
# Estimate Std. Error    z value     Pr(>|z|)
# (Intercept) 0.006403353  0.1052792 0.06082258 9.515005e-01
# mrs         0.613037418  0.1405462 4.36182285 1.289833e-05
# sva1        9.480318338  2.6183605 3.62070788 2.937981e-04
# sva2        7.591637312  2.2264799 3.40970402 6.503342e-04
# 
# $basic
# Estimate Std. Error    z value     Pr(>|z|)
# (Intercept)   3.75180146 1.38196189  2.7148371 6.630837e-03
# mrs           0.80156269 0.16367197  4.8973730 9.712636e-07
# age          -0.03479868 0.01943792 -1.7902470 7.341421e-02
# sex           0.43200824 0.26100760  1.6551558 9.789290e-02
# CD4T         -7.68616733 2.23797135 -3.4344351 5.937903e-04
# CD8T        -10.31984223 4.12088765 -2.5042765 1.227021e-02
# NK          -18.73136402 3.43328807 -5.4558090 4.875037e-08
# Bcell         3.13225615 5.49897685  0.5696071 5.689442e-01
# Mono          8.37504920 4.70030792  1.7818086 7.478045e-02
# sva1         11.81577539 2.95608679  3.9971003 6.412313e-05
# sva2         -2.20187565 2.99004346 -0.7364026 4.614857e-01
# 
# $plus_risk_factors
# Estimate Std. Error    z value     Pr(>|z|)
# (Intercept)   2.21173181 2.02036591  1.0947184 2.736400e-01
# mrs           0.50160823 0.19913039  2.5189939 1.176907e-02
# age          -0.01537427 0.02504326 -0.6139087 5.392757e-01
# sex          -0.15239478 0.31868461 -0.4781994 6.325083e-01
# CD4T         -8.99254642 2.70158058 -3.3286242 8.727609e-04
# CD8T        -10.59668142 5.21815964 -2.0307315 4.228223e-02
# NK          -17.76831765 4.46108190 -3.9829615 6.806177e-05
# Bcell         6.79272743 6.69901597  1.0139888 3.105881e-01
# Mono         14.61321480 5.79630244  2.5211270 1.169796e-02
# bmi          -0.04515969 0.03085969 -1.4633876 1.433613e-01
# diabetes      0.98154193 0.35830672  2.7393903 6.155326e-03
# smk_now       1.65000081 0.40277302  4.0966021 4.192586e-05
# h_lipid       0.56439315 0.28764233  1.9621352 4.974674e-02
# htn           0.21768844 0.30694423  0.7092117 4.781931e-01
# sva1          9.13186893 3.56093593  2.5644575 1.033372e-02
# sva2         -4.46463423 3.50135425 -1.2751164 2.022680e-01

combined_model_results
# $unadjusted
# Estimate Std. Error    z value    Pr(>|z|)
# (Intercept) 0.009970967  0.1045037 0.09541257 0.923987133
# mrs         0.509680130  0.1380247 3.69267335 0.000221909
# sva1        8.254259641  2.5840557 3.19430408 0.001401684
# sva2        6.085343031  2.1267347 2.86135498 0.004218344
# 
# $basic
# Estimate Std. Error   z value     Pr(>|z|)
# (Intercept)   3.69125204 1.41726610  2.604488 9.201174e-03
# mrs           0.55857307 0.17605058  3.172799 1.509771e-03
# age          -0.03988137 0.02067082 -1.929356 5.368663e-02
# sex           0.27675167 0.25537288  1.083716 2.784907e-01
# CD4T         -5.56492439 2.15979995 -2.576593 9.977950e-03
# CD8T         -9.71792476 4.07960845 -2.382073 1.721549e-02
# NK          -17.05726067 3.37801442 -5.049493 4.429838e-07
# Bcell         2.51574206 5.39683894  0.466151 6.411074e-01
# Mono          9.36670154 4.58366289  2.043497 4.100326e-02
# sva1          9.18784322 3.01277688  3.049626 2.291264e-03
# sva2         -3.23916694 2.97643265 -1.088272 2.764753e-01
# 
# $plus_risk_factors
# Estimate Std. Error    z value     Pr(>|z|)
# (Intercept)   1.78971312 2.05186519  0.8722372 3.830790e-01
# mrs           0.31102581 0.22211711  1.4002785 1.614299e-01
# age          -0.01453999 0.02617189 -0.5555575 5.785134e-01
# sex          -0.25467483 0.31685303 -0.8037633 4.215337e-01
# CD4T         -7.65368577 2.64971424 -2.8884948 3.870905e-03
# CD8T        -10.13494625 5.25407731 -1.9289679 5.373485e-02
# NK          -16.85823386 4.42260376 -3.8118346 1.379392e-04
# Bcell         6.57017678 6.63786023  0.9898034 3.222702e-01
# Mono         14.99588303 5.76683546  2.6003660 9.312437e-03
# bmi          -0.04220355 0.03079856 -1.3703089 1.705905e-01
# diabetes      1.08927434 0.35409883  3.0761874 2.096660e-03
# smk_now       1.72890760 0.39980514  4.3243756 1.529645e-05
# h_lipid       0.57095090 0.28598339  1.9964477 4.588521e-02
# htn           0.24213859 0.30427124  0.7957985 4.261492e-01
# sva1          6.98340411 3.66268123  1.9066372 5.656759e-02
# sva2         -4.95630510 3.55213839 -1.3953018 1.629248e-01

# ROC curves
metadata_roc <- na.omit(metadata)  # To ensure sample sizes are consistent across models
mrs_only_model <- test_mrs(metadata_roc, metadata_roc$mrs, 
                           "case ~ mrs", return_fit=T)
mrs_only_ROC <- roc(metadata_roc$case, predict(mrs_only_model, type="response"))
# frs_model <- glm(case ~ frs, data=metadata_roc, family="binomial")
# frs_ROC <- roc(metadata_roc$case, predict(frs_model, type="response"))
rf_model <- glm(as.formula(paste0(
  "case ~ age + sex + CD4T + CD8T + NK + Bcell + Mono + ",
  # "bmi + diabetes + smk_now + ldl + hdl + sbp")), 
  "bmi + diabetes + smk_now + h_lipid + htn + sva1 + sva2")),
data=metadata_roc, family="binomial")
rf_ROC <- roc(metadata_roc$case, predict(rf_model, type="response"))
# combined_model <- test_mrs(metadata_roc, metadata_roc$mrs, 
#                            "case ~ mrs + frs", return_fit=T)
combined_model <- glm(as.formula(paste0(
  "case ~ mrs + age + sex + CD4T + CD8T + NK + Bcell + Mono + ",
  # "bmi + diabetes + smk_now + ldl + hdl + sbp")), 
  "bmi + diabetes + smk_now + h_lipid + htn + sva1 + sva2")),
data=metadata_roc, family="binomial")
combined_ROC <- roc(metadata_roc$case, predict(combined_model, type="response"))
delong_test <- roc.test(rf_ROC, combined_ROC)
saveRDS(delong_test, "delong_test.rds")
delong_test
# DeLong's test for two correlated ROC curves
# data:  rf_ROC and combined_ROC
# Z = -1.5405, p-value = 0.1234
# alternative hypothesis: true difference in AUC is not equal to 0
# sample estimates:
# AUC of roc1 AUC of roc2 
#   0.8243908   0.8358068

png("roc_curves.png")
plot(mrs_only_ROC, col="blue")
# plot(frs_ROC, col="red", add=T)
plot(combined_ROC, col="green", add=T)
legend("bottomright", 
       legend=paste0(c("MRS: auc=", 
                       # "RF: auc=", 
                       "MRS+RF: auc="), 
                     round(as.numeric(c(mrs_only_ROC$auc, 
                                        # rf_ROC$auc, 
                                        combined_ROC$auc)), 2)), 
       col=c("blue", "red", "green"), lty=1, lwd=3, cex=1.3)
dev.off()

# Age, sex, and Framingham Risk Score interactions
test_mrs_with_SE <- function(nmd, subset, covar_string) {
  subset_data <- nmd[subset, ]
  test_res <- test_mrs(subset_data, subset_data$mrs, 
                       paste("case ~ mrs +", covar_string),
                       return_fit=T)
  model_coefs <- summary(test_res)$coef["mrs", c("Estimate", "Std. Error")]
  tibble(OR_per_SD=exp(model_coefs["Estimate"]),
         SE_lower=exp(model_coefs["Estimate"] - model_coefs["Std. Error"]),
         SE_higher=exp(model_coefs["Estimate"] + model_coefs["Std. Error"]),
         N=nrow(subset_data))
}

#AFS# without svas
sex_strat <- map_dfr(unique(metadata$sex), function(s) {
  tryCatch(test_mrs_with_SE(metadata, metadata$sex == s, 
                            "age + CD4T + CD8T + NK + Bcell + Mono"),
           error=function(e) data.frame(OR_per_SD=NA, SE_lower=NA, SE_higher=NA, N=NA))
}, .id="Sex") %>%
  mutate(Sex=unique(metadata$sex))  # Ensure sex label order is correct

metadata$age_group <- cut(metadata$age,
                          breaks=c(0, 55, 65, 75, 100),
                          labels=c("0-55", "56-65", "66-75", "76+"))
age_group_strat <- map_dfr(levels(metadata$age_group), function(g) {
  tryCatch(test_mrs_with_SE(metadata, metadata$age_group == g,
                            "age + sex + CD4T + CD8T + NK + Bcell + Mono"),  
           error=function(e) data.frame(OR_per_SD=NA, SE_lower=NA, SE_higher=NA, N=NA))
}, .id="Age group") %>%
  mutate(`Age group`=levels(metadata$age_group))

metadata_roc$rf_risk <- predict(rf_model, type="response")
metadata_roc$rf_risk_group <- cut(metadata_roc$rf_risk, 4, labels=paste0("Q", 1:4))
rf_risk_group_strat <- map_dfr(levels(metadata_roc$rf_risk_group), function(g) {
  tryCatch(test_mrs_with_SE(metadata_roc, metadata_roc$rf_risk_group == g,
                            "age + sex + CD4T + CD8T + NK + Bcell + Mono"),
           error=function(e) data.frame(OR_per_SD=NA, SE_lower=NA, SE_higher=NA, N=NA))
}, .id="Risk factor-based risk quartile") %>%
  mutate(`Risk factor-based risk quartile`=levels(metadata_roc$rf_risk_group))

save("sex_strat", "age_group_strat", "rf_risk_group_strat",
     file="stratified_results_nosvas.RData")


#AFS# with svas
sex_strat <- map_dfr(unique(metadata$sex), function(s) {
  tryCatch(test_mrs_with_SE(metadata, metadata$sex == s, 
                            "age + CD4T + CD8T + NK + Bcell + Mono + sva1 + sva2"),
           error=function(e) data.frame(OR_per_SD=NA, SE_lower=NA, SE_higher=NA, N=NA))
}, .id="Sex") %>%
  mutate(Sex=unique(metadata$sex))  # Ensure sex label order is correct

metadata$age_group <- cut(metadata$age,
                          breaks=c(0, 55, 65, 75, 100),
                          labels=c("0-55", "56-65", "66-75", "76+"))
age_group_strat <- map_dfr(levels(metadata$age_group), function(g) {
  tryCatch(test_mrs_with_SE(metadata, metadata$age_group == g,
                            "age + sex + CD4T + CD8T + NK + Bcell + Mono + sva1 + sva2"),  
           error=function(e) data.frame(OR_per_SD=NA, SE_lower=NA, SE_higher=NA, N=NA))
}, .id="Age group") %>%
  mutate(`Age group`=levels(metadata$age_group))

metadata_roc$rf_risk <- predict(rf_model, type="response")
metadata_roc$rf_risk_group <- cut(metadata_roc$rf_risk, 4, labels=paste0("Q", 1:4))
rf_risk_group_strat <- map_dfr(levels(metadata_roc$rf_risk_group), function(g) {
  tryCatch(test_mrs_with_SE(metadata_roc, metadata_roc$rf_risk_group == g,
                            "age + sex + CD4T + CD8T + NK + Bcell + Mono + sva1 + sva2"),
           error=function(e) data.frame(OR_per_SD=NA, SE_lower=NA, SE_higher=NA, N=NA))
}, .id="Risk factor-based risk quartile") %>%
  mutate(`Risk factor-based risk quartile`=levels(metadata_roc$rf_risk_group))

save("sex_strat", "age_group_strat", "rf_risk_group_strat",
     file="stratified_results.RData")


###AFS### sex_strat svas
# # A tibble: 2 x 5
#     Sex OR_per_SD SE_lower SE_higher     N
#   <dbl>     <dbl>    <dbl>     <dbl> <int>
# 1     0      1.73     1.37      2.18   190
# 2     1      3.18     2.49      4.08   201

###AFS### sex_strat no svas
# # A tibble: 2 x 5
#     Sex OR_per_SD SE_lower SE_higher     N
#   <dbl>     <dbl>    <dbl>     <dbl> <int>
# 1     0      1.33     1.11      1.59   190
# 2     1      1.90     1.59      2.28   201

###AFS### age_group_strat svas
# # A tibble: 4 x 5
# `Age group`   OR_per_SD SE_lower SE_higher    N
#   <chr>           <dbl>    <dbl>     <dbl> <int>
# 1 0-55             6.79     3.91     11.8   64
# 2 56-65            3.43     2.57     4.57   161
# 3 66-75            1.29     1.03     1.61   166
# 4 76+              NA       NA       NA     NA

###AFS### age_group_strat no svas
# # A tibble: 4 x 5
#   `Age group` OR_per_SD SE_lower SE_higher     N
#   <chr>           <dbl>    <dbl>     <dbl> <int>
# 1 0-55             4.59    2.88       7.32    64
# 2 56-65            1.97    1.61       2.40   161
# 3 66-75            1.03    0.862      1.24   166
# 4 76+             NA      NA         NA       NA

###AFS### rf_risk_group_strat svas
# # A tibble: 4 x 5
# `Risk factor-based risk quartile`   OR_per_SD SE_lower SE_higher     N
# <chr>                                   <dbl>    <dbl>     <dbl> <int>
# 1 Q1                                     4.97    2.88       8.57   119
# 2 Q2                                     1.12    0.831      1.51    90
# 3 Q3                                     2.92    1.75       4.89    60
# 4 Q4                                     1.36    0.639      2.91    55

###AFS### rf_risk_group_strat no svas
# # A tibble: 4 x 5
# `Risk factor-based risk quartile`   OR_per_SD SE_lower SE_higher     N
# <chr>                                   <dbl>    <dbl>     <dbl> <int>
# 1 Q1                                    1.99     1.39      2.86    119
# 2 Q2                                    1.17     0.922     1.49     90
# 3 Q3                                    2.10     1.47      3.01     60
# 4 Q4                                    0.539    0.345     0.841    55