library(tidyverse)
library(survival)
betas_lbc <- readRDS("../int/betas.qc.norm.filt_lbc.rds")
metaData <- readRDS("../int/metaData.rds")
estCellCounts_lbc <- readRDS("../int/est_cell_counts_lbc.rds")
load("../int/CPACOR_lbc.RData")
CP_PCs_lbc <- CP_PCs
load("../int/PCA.fit_lbc.RData")
PCs_lbc <- PCs
nonMethData <- Reduce(function(x,y) inner_join(x, y, by="sampleKey"), 
                      list(metaData, estCellCounts_lbc, CP_PCs_lbc, PCs_lbc)) %>%
  filter(is.na(wave) | wave==1) %>%
  distinct(subjID, .keep_all=T) 

nmd_lbc21 <- filter(nonMethData, study=="lbc21", wave==1, sampleKey %in% colnames(betas_lbc))
betas_lbc21 <- betas_lbc[,nmd_lbc21$sampleKey]
nmd_lbc36 <- filter(nonMethData, study=="lbc36", wave==1, sampleKey %in% colnames(betas_lbc))
betas_lbc36 <- betas_lbc[,nmd_lbc36$sampleKey]

load("../cache/mrs_models/train-mcmm-full_8b2c597d335dbe4b109b61d82046fefa.RData")
calc_mrs_glmnet <- function(coefs, meth) {
  if (length(coefs)==0 || all(coefs==0)) {
    message("No non-zero coefficients")
    return(rep(NA, ncol(meth)))
  } else {
    nonZeroCoefs <- coefs[coefs!=0]
    usefulCoefs <- nonZeroCoefs[names(nonZeroCoefs) %in% rownames(meth)]
    as.vector(t(meth[names(usefulCoefs),,drop=F]) %*% usefulCoefs)
  }
} 
test_mrs <- function(meta, mrs, subset) {
  meta$mrs <- scale(mrs)
  meta$cvd_surv <- Surv(time=meta$time, event=meta$event)  # Outcome
  mrs.test <- coxph(cvd_surv~mrs, data=meta, subset=subset) # Test the MRS using a Cox model
  mrs.res <- summary(mrs.test)$coef["mrs",c("exp(coef)","z")]
  tibble(HR_per_SD=mrs.res["exp(coef)"], p=2*pnorm(-abs(mrs.res["z"])))
}
calc_test_mrs <- function(mrs_fit, betas, meta, subset=NULL) {
  coefs <- mrs_fit$cpgCoefs
  mrs <- calc_mrs_glmnet(coefs, betas)
  test_mrs(meta, mrs, subset)
}

mrs21 <- calc_mrs_glmnet(mrs.fit_WhiFhs_full$cpgCoefs, betas_lbc21)
mrs36 <- calc_mrs_glmnet(mrs.fit_WhiFhs_full$cpgCoefs, betas_lbc36)