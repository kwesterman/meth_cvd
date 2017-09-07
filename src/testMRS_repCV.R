suppressMessages(silent <- lapply(c("tidyverse","survival","gtools","glmnet","minfi","caret"), 
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
complete_cases <- complete.cases(covars_mat)
design_mat <- cbind(covars_mat, t(Mvals))[complete_cases,]
rm(Mvals)

# Additional needed variables
shareids <- nonMethData$shareid[complete_cases]
cvd_surv <- Surv(time=nonMethData$timeToEvent, event=nonMethData$event)[complete_cases]  # Outcome
frs <- nonMethData$frs[complete_cases]
pastEvent <- nonMethData$pastEvent[complete_cases]

trainTestMRS <- function(cvd_surv, design_mat, frs, pastEvent, shareids, trainset_idx) {
  train_ids <- shareids[trainset_idx]
  mrs.fit <- glmnet(design_mat[trainset_idx,], cvd_surv[trainset_idx], family="cox", alpha=0.5,
                    penalty.factor=ifelse(grepl("^cg", colnames(design_mat)), 1, 0))
  coefs <- mrs.fit$beta[,ncol(mrs.fit$beta)]
  coefs_meth <- coefs[coefs!=0 & grepl("cg", names(coefs))]
  mrs_calc <- as.vector(design_mat[,names(coefs_meth)] %*% coefs_meth)
  testset <- setdiff(1:nrow(design_mat), trainset_idx)
  mrs.test <- coxph(cvd_surv~mrs, data=data.frame(cvd_surv=cvd_surv, mrs=mrs_calc), subset=testset)
  mrs_pastEvent.test <- coxph(cvd_surv~mrs+pastEvent, 
                        data=data.frame(cvd_surv=cvd_surv, mrs=mrs_calc, pastEvent=pastEvent), 
                        subset=testset)
  mrs_frs.test <- coxph(cvd_surv~mrs+frs, 
                        data=data.frame(cvd_surv=cvd_surv, mrs=mrs_calc, frs=frs), 
                        subset=testset)
  list(components=coefs_meth, 
       res=summary(mrs.test)$coef, 
       res_withPastEvent=summary(mrs_pastEvent.test)$coef,
       res_withFRS=summary(mrs_frs.test)$coef, 
       train_ids=train_ids)
}

set.seed(1)
splits <- createDataPartition(factor(cvd_surv[,"status"]), times=20, p=0.6)  # Set of train/test splits w/ 
allRes <- lapply(1:20, function(idx) {
  print(paste("Iteration", idx))
  trainTestMRS(cvd_surv, design_mat, frs, pastEvent, shareids, splits[[idx]])
})

save("allRes", file="../int/testMRS_repCV_res_1CPA_noPastEvents.RData")


