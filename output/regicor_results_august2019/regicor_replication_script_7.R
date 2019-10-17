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

setwd("/projects/regicor/METIAM/EPIC/misc/kenny_replication/sept2019/")

#AFS##Load betas
load("../betas.RData")

#AFS##Load metadata
load("../2019feb/metadata.RData")


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

# Setup and MRS calculation
stopifnot(all(metadata$subject_id == colnames(betas)))
ssl_coefs <- readRDS("ssl_mrs_coefs.rds")  # List of CpG weight vectors per study
study_weights <- readRDS("csl_weights.rds")  # Named vector of study weights
metadata$mrs <- assemble_csl(metadata, betas, ssl_coefs, study_weights)$mrs
combined_model_coefs <- readRDS("combined_mrs_coefs.rds")  # Named vector of coefficient weights
metadata$mrs_combined <- calc_mrs(combined_model_coefs, betas)
combat_model_coefs <- readRDS("combat_mrs_coefs.rds")
metadata$mrs_combat <- calc_mrs(combat_model_coefs, betas)

# Initial evaluation of MRS performance
# (Note: assumes case/control variable is called "case" -- alter as needed)
# Series of models to test using test_mrs as above:
model_strings <- list(
  unadjusted="case ~ mrs + sva1 + sva2",
  basic="case ~ mrs + age + sex + CD4T + CD8T + NK + Bcell + Mono + sva1 + sva2",
  plus_risk_factors=paste0(
    "case ~ mrs + age + sex + CD4T + CD8T + NK + Bcell + Mono + ",
    "bmi + diabetes + smk_now + h_lipid + htn + sva1 + sva2")
)
csl_model_results <- lapply(model_strings, function(ms) {
  test_mrs(metadata, metadata$mrs, ms)
})
combined_model_results <- lapply(
  model_strings, function(ms) test_mrs(metadata, metadata$mrs_combined, ms))
combat_model_results <- lapply(
  model_strings, function(ms) test_mrs(metadata, metadata$mrs_combat, ms))
save("csl_model_results", "combined_model_results", "combat_model_results",
     file="main_mrs_results.RData")

csl_model_results
# $unadjusted
#                Estimate Std. Error    z value     Pr(>|z|)
# (Intercept) 0.005487968  0.1054814 0.05202781 9.585065e-01
# mrs         0.602612601  0.1334681 4.51503191 6.330723e-06
# sva1        8.713167117  2.4929215 3.49516304 4.737722e-04
# sva2        7.262676372  2.2073383 3.29024168 1.001014e-03
# 
# $basic
#                 Estimate Std. Error    z value     Pr(>|z|)
# (Intercept)   3.48788777 1.36670614  2.5520393 1.070944e-02
# mrs           0.75952969 0.15466761  4.9107223 9.074152e-07
# age          -0.03202376 0.01930744 -1.6586225 9.719189e-02
# sex           0.46504858 0.26187901  1.7758147 7.576347e-02
# CD4T         -7.39616379 2.22995543 -3.3167317 9.107704e-04
# CD8T        -10.03645609 4.12025397 -2.4358829 1.485549e-02
# NK          -18.84402404 3.43704643 -5.4826213 4.190694e-08
# Bcell         3.48189766 5.54017303  0.6284818 5.296884e-01
# Mono          8.48098209 4.70226546  1.8035949 7.129483e-02
# sva1         10.46537739 2.78763196  3.7542178 1.738837e-04
# sva2         -2.59764189 2.97504703 -0.8731431 3.825851e-01
# 
# $plus_risk_factors
#                 Estimate Std. Error    z value     Pr(>|z|)
# (Intercept)   2.08648145 2.00740114  1.0393944 2.986214e-01
# mrs           0.47789590 0.19005643  2.5144947 1.192031e-02
# age          -0.01369318 0.02489834 -0.5499637 5.823443e-01
# sex          -0.12950802 0.31930424 -0.4055944 6.850406e-01
# CD4T         -8.85959249 2.69934703 -3.2821243 1.030282e-03
# CD8T        -10.53567236 5.21671708 -2.0195982 4.342508e-02
# NK          -17.85143272 4.46362256 -3.9993150 6.352608e-05
# Bcell         6.94321483 6.71346013  1.0342230 3.010319e-01
# Mono         14.64903358 5.79463820  2.5280325 1.147037e-02
# bmi          -0.04557140 0.03091865 -1.4739130 1.405050e-01
# diabetes      0.97707776 0.35982403  2.7154322 6.618933e-03
# smk_now       1.64755817 0.40213908  4.0969860 4.185643e-05
# h_lipid       0.55725351 0.28762363  1.9374399 5.269159e-02
# htn           0.21911391 0.30682516  0.7141328 4.751451e-01
# sva1          8.31208759 3.36538070  2.4698803 1.351583e-02
# sva2         -4.74617440 3.48137237 -1.3633056 1.727862e-01
combined_model_results
# $unadjusted
#               Estimate Std. Error   z value     Pr(>|z|)
# (Intercept) 0.01328912  0.1061950 0.1251388 9.004136e-01
# mrs         0.61906149  0.1266986 4.8860941 1.028560e-06
# sva1        6.86453224  2.2823114 3.0077106 2.632237e-03
# sva2        5.47709188  2.1477746 2.5501241 1.076846e-02
# 
# $basic
#                 Estimate Std. Error    z value     Pr(>|z|)
# (Intercept)   4.42930763 1.42270435  3.1133015 1.850069e-03
# mrs           0.75273456 0.15403174  4.8868795 1.024467e-06
# age          -0.05082421 0.02050504 -2.4786205 1.318916e-02
# sex           0.43071893 0.26234003  1.6418346 1.006243e-01
# CD4T         -6.58547970 2.20777170 -2.9828626 2.855661e-03
# CD8T         -8.97751339 4.17476661 -2.1504228 3.152178e-02
# NK          -17.69195083 3.45809414 -5.1160987 3.119199e-07
# Bcell         3.76432125 5.53790262  0.6797377 4.966706e-01
# Mono          8.70601886 4.69926423  1.8526345 6.393477e-02
# sva1          7.95063650 2.54063615  3.1293881 1.751708e-03
# sva2         -3.65245076 2.94577834 -1.2398933 2.150149e-01
# 
# $plus_risk_factors
#                 Estimate Std. Error    z value     Pr(>|z|)
# (Intercept)   2.91314391 2.10503580  1.3838928 1.663913e-01
# mrs           0.50849083 0.18623768  2.7303327 6.327044e-03
# age          -0.02694410 0.02626947 -1.0256813 3.050418e-01
# sex          -0.16399304 0.31972188 -0.5129240 6.080045e-01
# CD4T         -8.44174305 2.67719516 -3.1532042 1.614888e-03
# CD8T         -9.71621382 5.29231322 -1.8359106 6.637088e-02
# NK          -17.62733611 4.52159084 -3.8984810 9.679798e-05
# Bcell         7.40217526 6.75048395  1.0965399 2.728426e-01
# Mono         14.43613473 5.81704763  2.4816944 1.307594e-02
# bmi          -0.04798892 0.03138993 -1.5287997 1.263141e-01
# diabetes      1.02663140 0.35825021  2.8656826 4.161112e-03
# smk_now       1.59425241 0.40258931  3.9599969 7.495074e-05
# h_lipid       0.52057498 0.28883730  1.8023122 7.149630e-02
# htn           0.20340539 0.30729479  0.6619227 5.080208e-01
# sva1          7.07575084 3.06224096  2.3106447 2.085249e-02
# sva2         -5.26225574 3.46931950 -1.5167977 1.293178e-01
combat_model_results
# $unadjusted
#               Estimate Std. Error   z value     Pr(>|z|)
# (Intercept) 0.01176699  0.1055690 0.1114626 9.112496e-01
# mrs         0.58152104  0.1296464 4.4854381 7.276426e-06
# sva1        7.45835703  2.3559450 3.1657602 1.546782e-03
# sva2        5.26416553  2.1332703 2.4676506 1.360030e-02
# 
# $basic
#                 Estimate Std. Error    z value     Pr(>|z|)
# (Intercept)   4.51008862 1.42689881  3.1607628 1.573566e-03
# mrs           0.76779542 0.15750570  4.8747151 1.089656e-06
# age          -0.04942207 0.02033643 -2.4302232 1.508953e-02
# sex           0.46919252 0.26244027  1.7878069 7.380717e-02
# CD4T         -7.44825043 2.23478994 -3.3328638 8.595702e-04
# CD8T        -10.31299978 4.15817768 -2.4801729 1.313187e-02
# NK          -17.81053094 3.45232836 -5.1589910 2.482843e-07
# Bcell         3.28775805 5.54137873  0.5933105 5.529734e-01
# Mono          9.04491287 4.71919850  1.9166206 5.528615e-02
# sva1          9.22320469 2.65698768  3.4713013 5.179426e-04
# sva2         -4.71957873 2.93362680 -1.6087863 1.076631e-01
# 
# $plus_risk_factors
#                 Estimate Std. Error    z value     Pr(>|z|)
# (Intercept)   3.06162180 2.09740926  1.4597160 1.443681e-01
# mrs           0.56579210 0.18669202  3.0306175 2.440542e-03
# age          -0.02822745 0.02603090 -1.0843823 2.781953e-01
# sex          -0.13273603 0.32081745 -0.4137432 6.790622e-01
# CD4T         -9.08154552 2.72686104 -3.3304028 8.672045e-04
# CD8T        -10.60848263 5.28311094 -2.0079992 4.464337e-02
# NK          -17.60774644 4.52847711 -3.8882269 1.009792e-04
# Bcell         7.16756391 6.80208649  1.0537302 2.920065e-01
# Mono         14.92067029 5.83741909  2.5560389 1.058713e-02
# bmi          -0.04783578 0.03154911 -1.5162324 1.294606e-01
# diabetes      0.99373032 0.35913351  2.7670220 5.657094e-03
# smk_now       1.60882914 0.40427798  3.9795121 6.905683e-05
# h_lipid       0.52659432 0.28910909  1.8214381 6.854028e-02
# htn           0.22462885 0.30726326  0.7310632 4.647406e-01
# sva1          8.22555665 3.19400291  2.5753128 1.001495e-02
# sva2         -5.99414106 3.45900168 -1.7329107 8.311154e-02

# ROC curves
metadata_roc <- na.omit(metadata)  # To ensure sample sizes are consistent across models
mrs_only_model <- test_mrs(metadata_roc, metadata_roc$mrs, 
                           "case ~ mrs", return_fit=T)
mrs_only_ROC <- roc(metadata_roc$case, predict(mrs_only_model, type="response"))
rf_model <- glm(as.formula(paste0(
  "case ~ age + sex + CD4T + CD8T + NK + Bcell + Mono + ",
  "bmi + diabetes + smk_now + h_lipid + htn + sva1 + sva2")),
data=metadata_roc, family="binomial")
rf_ROC <- roc(metadata_roc$case, predict(rf_model, type="response"))
combined_model <- glm(as.formula(paste0(
  "case ~ mrs + age + sex + CD4T + CD8T + NK + Bcell + Mono + ",
  "bmi + diabetes + smk_now + h_lipid + htn + sva1 + sva2")),
data=metadata_roc, family="binomial")
combined_ROC <- roc(metadata_roc$case, predict(combined_model, type="response"))
delong_test <- roc.test(rf_ROC, combined_ROC)
saveRDS(delong_test, "delong_test.rds")

delong_test
# DeLong's test for two correlated ROC curves
# 
# data:  rf_ROC and combined_ROC
# Z = -1.5815, p-value = 0.1138
# alternative hypothesis: true difference in AUC is not equal to 0
# sample estimates:
# AUC of roc1 AUC of roc2 
# 0.8243908   0.8363185 

png("roc_curves.png")
plot(mrs_only_ROC, col="blue")
plot(combined_ROC, col="green", add=T)
legend("bottomright", 
       legend=paste0(c("MRS: auc=", 
                       "MRS+RF: auc="), 
                     round(as.numeric(c(mrs_only_ROC$auc, 
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

metadata_roc$rf_risk <- predict(rf_model, type="response")
metadata_roc$rf_risk_group <- cut(metadata_roc$rf_risk, 4, labels=paste0("Q", 1:4))
rf_risk_group_strat_csl <- map_dfr(levels(metadata_roc$rf_risk_group), function(g) {
  tryCatch(test_mrs_with_SE(metadata_roc, metadata_roc$rf_risk_group == g,
                            "age + sex + CD4T + CD8T + NK + Bcell + Mono + sva1 + sva2"),
           error=function(e) data.frame(OR_per_SD=NA, SE_lower=NA, SE_higher=NA, N=NA))
}, .id="Risk factor-based risk quartile") %>%
  mutate(`Risk factor-based risk quartile`=levels(metadata_roc$rf_risk_group))
rf_risk_group_strat_combined <- map_dfr(levels(metadata_roc$rf_risk_group), function(g) {
  tryCatch(test_mrs_with_SE(mutate(metadata_roc, mrs=mrs_combined),
                            metadata_roc$rf_risk_group == g,
                            "age + sex + CD4T + CD8T + NK + Bcell + Mono + sva1 + sva2"),
           error=function(e) data.frame(OR_per_SD=NA, SE_lower=NA, SE_higher=NA, N=NA))
}, .id="Risk factor-based risk quartile") %>%
  mutate(`Risk factor-based risk quartile`=levels(metadata_roc$rf_risk_group))
rf_risk_group_strat_combat <- map_dfr(levels(metadata_roc$rf_risk_group), function(g) {
  tryCatch(test_mrs_with_SE(mutate(metadata_roc, mrs=mrs_combat),
                            metadata_roc$rf_risk_group == g,
                            "age + sex + CD4T + CD8T + NK + Bcell + Mono + sva1 + sva2"),
           error=function(e) data.frame(OR_per_SD=NA, SE_lower=NA, SE_higher=NA, N=NA))
}, .id="Risk factor-based risk quartile") %>%
  mutate(`Risk factor-based risk quartile`=levels(metadata_roc$rf_risk_group))

save(
  # "sex_strat", "age_group_strat", 
  "rf_risk_group_strat_csl",
  "rf_risk_group_strat_combined", "rf_risk_group_strat_combat",
  file="stratified_results.RData"
)

## OLD -- UNNEEDED AT THE MOMENT

#AFS# without svas

# sex_strat <- map_dfr(unique(metadata$sex), function(s) {
#   tryCatch(test_mrs_with_SE(metadata, metadata$sex == s, 
#                             "age + CD4T + CD8T + NK + Bcell + Mono"),
#            error=function(e) data.frame(OR_per_SD=NA, SE_lower=NA, SE_higher=NA, N=NA))
# }, .id="Sex") %>%
#   mutate(Sex=unique(metadata$sex))  # Ensure sex label order is correct
# 
# metadata$age_group <- cut(metadata$age,
#                           breaks=c(0, 55, 65, 75, 100),
#                           labels=c("0-55", "56-65", "66-75", "76+"))
# age_group_strat <- map_dfr(levels(metadata$age_group), function(g) {
#   tryCatch(test_mrs_with_SE(metadata, metadata$age_group == g,
#                             "age + sex + CD4T + CD8T + NK + Bcell + Mono"),  
#            error=function(e) data.frame(OR_per_SD=NA, SE_lower=NA, SE_higher=NA, N=NA))
# }, .id="Age group") %>%
#   mutate(`Age group`=levels(metadata$age_group))
# 
# metadata_roc$rf_risk <- predict(rf_model, type="response")
# metadata_roc$rf_risk_group <- cut(metadata_roc$rf_risk, 4, labels=paste0("Q", 1:4))
# rf_risk_group_strat <- map_dfr(levels(metadata_roc$rf_risk_group), function(g) {
#   tryCatch(test_mrs_with_SE(metadata_roc, metadata_roc$rf_risk_group == g,
#                             "age + sex + CD4T + CD8T + NK + Bcell + Mono"),
#            error=function(e) data.frame(OR_per_SD=NA, SE_lower=NA, SE_higher=NA, N=NA))
# }, .id="Risk factor-based risk quartile") %>%
#   mutate(`Risk factor-based risk quartile`=levels(metadata_roc$rf_risk_group))
# 
# save("sex_strat", "age_group_strat", "rf_risk_group_strat",
#      file="stratified_results_nosvas.RData")


#AFS# with svas

# sex_strat <- map_dfr(unique(metadata$sex), function(s) {
#   tryCatch(test_mrs_with_SE(metadata, metadata$sex == s, 
#                             "age + CD4T + CD8T + NK + Bcell + Mono + sva1 + sva2"),
#            error=function(e) data.frame(OR_per_SD=NA, SE_lower=NA, SE_higher=NA, N=NA))
# }, .id="Sex") %>%
#   mutate(Sex=unique(metadata$sex))  # Ensure sex label order is correct
# 
# metadata$age_group <- cut(metadata$age,
#                           breaks=c(0, 55, 65, 75, 100),
#                           labels=c("0-55", "56-65", "66-75", "76+"))
# age_group_strat <- map_dfr(levels(metadata$age_group), function(g) {
#   tryCatch(test_mrs_with_SE(metadata, metadata$age_group == g,
#                             "age + sex + CD4T + CD8T + NK + Bcell + Mono + sva1 + sva2"),  
#            error=function(e) data.frame(OR_per_SD=NA, SE_lower=NA, SE_higher=NA, N=NA))
# }, .id="Age group") %>%
#   mutate(`Age group`=levels(metadata$age_group))
