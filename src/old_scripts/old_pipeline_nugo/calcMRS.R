suppressMessages(silent <- lapply(c("methods","tidyverse","survival","glmnet","caret"), 
                                  library, character.only=T))

# Load methylation data
load("../int/Mvals.RData")  # Loads Mvals object
print("Methylation data loaded.")

# Load non-methylation (event + covariate) data
load("../int/nonMethData.RData")
nonMethData <- nonMethData %>% replace_na(list(pastEvent=FALSE)) %>% filter(pastEvent==FALSE)

# Make sure samples and their ordering are identical in methylation and metadata
Mvals <- Mvals[,match(nonMethData$sampleKey, colnames(Mvals))]
print(paste("Dimensions of event/covariate matrix:", paste(dim(nonMethData), collapse=" x ")))
print(paste("Dimensions of M-value matrix:", paste(dim(Mvals), collapse=" x ")))
stopifnot(all(colnames(Mvals)==nonMethData$sampleKey))  # Ensure identical number and order of samples for methylation and covariate data

# Variables for regressions (including methylation M-values)
Mval_variances <- apply(Mvals, 1, var)
Mvals <- Mvals[Mval_variances>0.5,]

covars_formula <- paste0("~sex+age+smoking_now+bmi+CD4T+NK+Bcell+Mono+Gran+PC1_cp")
covars_frame <- model.frame(as.formula(covars_formula), nonMethData, na.action=na.pass)  # To allow NAs to pass through
covars_mat <- model.matrix(as.formula(covars_formula), covars_frame)[,-1]  # Converts factors to numeric
complete_cases <- complete.cases(covars_mat)  # Can't have missing values in glmnet input
design_mat <- cbind(covars_mat, t(Mvals))[complete_cases,]
cvd_surv <- Surv(time=nonMethData$timeToEvent, event=nonMethData$event)[complete_cases]  # Outcome

# Perform elastic net regression to generate MRS coefficients
print("Generating full MRS...")
mrs.fit <- glmnet(design_mat, cvd_surv, family="cox", alpha=0.5,
                  penalty.factor=ifelse(grepl("^cg", colnames(design_mat)), 1, 0))
coefs <- mrs.fit$beta[,ncol(mrs.fit$beta)]
mrs_coefs <- coefs[coefs!=0 & grepl("cg", names(coefs))]

# Repeat elastic net regression on a subset of 50% of samples
print("Generating half MRS")
set.seed(10)
train_idx <- createFolds(factor(cvd_surv[,"status"]), 2)[[1]]
train_ids <- nonMethData$shareid[complete_cases][train_idx]
mrs.fit <- glmnet(design_mat[train_idx,], cvd_surv[train_idx], family="cox", alpha=0.5,
                  penalty.factor=ifelse(grepl("^cg", colnames(design_mat)), 1, 0))
coefs <- mrs.fit$beta[,ncol(mrs.fit$beta)]
mrs_coefs_half <- coefs[coefs!=0 & grepl("cg", names(coefs))]

# Calculate MRS for all subjects
mrs_calc <- as.vector(t(Mvals)[,names(mrs_coefs)] %*% mrs_coefs)
mrs_calc_half <- as.vector(t(Mvals)[,names(mrs_coefs_half)] %*% mrs_coefs_half)

# Return final dataset
mrsData <- nonMethData %>%
  dplyr::select(shareid, CD8T:Gran, one_of(paste0("PC",1:20)),
                event, timeToEvent, pastEvent, 
                sex, age, bmi, smoking_now, CVD_med, HT_med, T2D_med,
                SBP, TOT_CHOL, HDL_CHOL, TRIG, GLUCOSE, frs) %>%
  dplyr::rename(sysBP=SBP, chol=TOT_CHOL, hdl=HDL_CHOL, tg=TRIG, glu=GLUCOSE) %>%
  dplyr::mutate(mrs=mrs_calc,
                mrs_half=mrs_calc_half,
                inTrainSet=shareid %in% train_ids)

save("mrsData", "mrs_coefs", file="../int/mrsData_1CPA_noPastEvents.RData")
