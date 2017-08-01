suppressMessages(silent <- lapply(c("methods","tidyverse","survival","glmnet","caret"), 
                                  library, character.only=T))

# Load methylation data
load("../int/Mvals.RData")  # Loads Mvals object
print("Methylation data loaded.")

# Load non-methylation (event + covariate) data
load("../int/nonMethData.RData")

# Sanity checks
final_id_set <- intersect(nonMethData$sampleKey, colnames(Mvals))
nonMethData <- nonMethData[match(final_id_set, nonMethData$sampleKey),]
Mvals <- Mvals[,final_id_set]
print(paste("Dimensions of event/covariate matrix:", paste(dim(nonMethData), collapse=" x ")))
print(paste("Dimensions of M-value matrix:", paste(dim(Mvals), collapse=" x ")))
stopifnot(ncol(Mvals)==nrow(nonMethData),  # Ensure same # of samples in methylation and covariate data
          all(colnames(Mvals)==nonMethData$sampleKey))  # Ensure identical order of samples for methylation and covariate data

# Variables for regressions (including methylation M-values)
covars_formula <- paste0("~sex+age+smoking_now+bmi+CD8T+CD4T+NK+Bcell+Mono+Gran+", 
                         paste0("PC",1:20,collapse="+"))
covars_mat <- data.matrix(model.frame(as.formula(covars_formula), nonMethData, na.action=na.pass))
complete_cases <- complete.cases(covars_mat)  # Can't have missing values in glmnet input
design_mat <- cbind(covars_mat, t(Mvals))[complete_cases,]
cvd_surv <- Surv(time=nonMethData$timeToEvent, event=nonMethData$event)[complete_cases]  # Outcome

# Perform elastic net regression to generate MRS coefficients
mrs.fit <- glmnet(design_mat, cvd_surv, family="cox", alpha=0.5)
coefs <- mrs.fit$beta[,ncol(mrs.fit$beta)]
coefs_meth_full <- coefs[coefs!=0 & grepl("cg", names(coefs))]

# Repeat elastic net regression on a subset of 50% of samples
set.seed(50)
train_idx <- createFolds(factor(cvd_surv[,"status"]), 2)[[1]]
train_ids <- nonMethData$shareid[complete_cases][train_idx]
mrs.fit <- glmnet(design_mat[train_idx,], cvd_surv[train_idx], family="cox", alpha=0.5)
coefs <- mrs.fit$beta[,ncol(mrs.fit$beta)]
coefs_meth_half <- coefs[coefs!=0 & grepl("cg", names(coefs))]

# Calculate MRS for all subjects
mrs_calc_full <- as.vector(t(Mvals)[,names(coefs_meth_full)] %*% coefs_meth_full)
mrs_calc_half <- as.vector(t(Mvals)[,names(coefs_meth_half)] %*% coefs_meth_half)

# Return final dataset
mrsData <- nonMethData %>%
  dplyr::select(shareid, CD8T:Gran, one_of(paste0("PC",1:20)),
                sex, age, bmi, smoking_now, CVD_med, HT_med, T2D_med,
                SBP, TOT_CHOL, HDL_CHOL, TRIG, GLUCOSE) %>%
  dplyr::rename(sysBP=SBP, chol=TOT_CHOL, hdl=HDL_CHOL, tg=TRIG, glu=GLUCOSE) %>%
  dplyr::mutate(mrs_full=mrs_calc_full,
                mrs_half=mrs_calc_half,
                inTrainSet=shareid %in% train_ids)

save("mrsData", file="../int/mrsData.RData")
