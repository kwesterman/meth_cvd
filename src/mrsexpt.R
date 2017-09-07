suppressMessages(silent <- lapply(c("tidyverse","survival","gtools","glmnet","minfi","caret"), 
                                  library, character.only=T))

# Load methylation data
load("../int/Mvals.RData")  # Loads Mvals object
print("Methylation data loaded.")

# Load non-methylation (event + covariate) data
load("../int/nonMethData.RData")

# Make sure samples and their ordering are identical in methylation and metadata
Mvals <- Mvals[,match(nonMethData$sampleKey, colnames(Mvals))]
print(paste("Dimensions of event/covariate matrix:", paste(dim(nonMethData), collapse=" x ")))
print(paste("Dimensions of M-value matrix:", paste(dim(Mvals), collapse=" x ")))
stopifnot(all(colnames(Mvals)==nonMethData$sampleKey))  # Ensure identical number and order of samples for methylation and covariate data

# Covariates for regressions (including methylation M-values)
models <- c(
  paste0("~sex+age+smoking_now+bmi+CD4T+NK+Bcell+Mono+Gran+", paste0("PC", 1:3, "_cp", collapse="+")),
  paste0("~sex+age+smoking_now+bmi+CD4T+NK+Bcell+Mono+Gran+", paste0("PC", 1:5, "_cp", collapse="+"))
)

trainTestMRS <- function(mod) {
  covars_frame <- model.frame(as.formula(mod), nonMethData, na.action=na.pass)  # To allow NAs to pass through
  covars_mat <- model.matrix(as.formula(mod), covars_frame)[,-1]  # Converts factors to numeric
  complete_cases <- complete.cases(covars_mat)
  design_mat <- cbind(covars_mat, t(Mvals))[complete_cases,]
  cvd_surv <- Surv(time=nonMethData$timeToEvent, event=nonMethData$event)[complete_cases]  # Outcome
  
  set.seed(25)
  trainset_idx <- createDataPartition(factor(cvd_surv[,"status"]), p=0.6)[[1]]
  
  mrs.fit <- glmnet(design_mat[trainset_idx,], cvd_surv[trainset_idx], family="cox", alpha=0.5,
                    penalty.factor=ifelse(grepl("^cg", colnames(design_mat)), 1, 0))
  coefs <- mrs.fit$beta[,ncol(mrs.fit$beta)]
  coefs_meth <- coefs[coefs!=0 & grepl("cg", names(coefs))]
  mrs_calc <- as.vector(design_mat[,names(coefs_meth)] %*% coefs_meth)
  testset <- setdiff(1:nrow(design_mat), trainset_idx)
  mrs.test <- coxph(cvd_surv~mrs, data=cbind(cvd_surv=cvd_surv, mrs=mrs_calc, 
                                             nonMethData[complete_cases,]), subset=testset)
  print(summary(mrs.test)$coef)
  mrs.test <- coxph(cvd_surv~mrs+frs, data=cbind(cvd_surv=cvd_surv, mrs=mrs_calc, 
                                             nonMethData[complete_cases,]), subset=testset)
  print(summary(mrs.test)$coef)
  mrs.test <- coxph(cvd_surv~mrs+pastEvent, data=cbind(cvd_surv=cvd_surv, mrs=mrs_calc, 
                                             nonMethData[complete_cases,]), subset=testset)
  print(summary(mrs.test)$coef)
  mrs.test <- coxph(cvd_surv~mrs+frs+pastEvent, data=cbind(cvd_surv=cvd_surv, mrs=mrs_calc, 
                                             nonMethData[complete_cases,]), subset=testset)
  print(summary(mrs.test)$coef)
  list(components=coefs_meth, res=summary(mrs.test)$coef)
}


allRes <- lapply(models, function(mod) {
  print(mod)
  trainTestMRS(mod)
})

save("allRes", file="../int/addPCs.RData")


