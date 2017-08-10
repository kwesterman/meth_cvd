suppressMessages(silent <- lapply(c("tidyverse","survival","gtools","glmnet","minfi","caret"), 
                                  library, character.only=T))
source("helpers.R")

### CREATE MRS

# Load methylation and sample data
load("../int/Mvals.RData")  # Loads Mvals object
print("Methylation data loaded.")

# Load non-methylation (event + covariate) data and sample similarity matrix
load("../int/nonMethData.RData")

# Sanity checks
final_id_set <- intersect(nonMethData$sampleKey, colnames(Mvals))
nonMethData <- nonMethData[match(final_id_set, nonMethData$sampleKey),]
Mvals <- Mvals[,final_id_set]
print(paste("Dimensions of event/covariate matrix:", paste(dim(nonMethData), collapse=" x ")))
print(paste("Dimensions of M-value matrix:", paste(dim(Mvals), collapse=" x ")))
stopifnot(ncol(Mvals)==nrow(nonMethData),  # Ensure same # of samples in methylation and covariate data
          all(colnames(Mvals)==nonMethData$sampleKey))  # Ensure identical order of samples for methylation and covariate data

# Covariates for regressions (including methylation M-values)
Bvals <- ilogit2(Mvals)
covars_formula <- paste0("~sex+age+smoking_now+bmi+CD8T+CD4T+NK+Bcell+Mono+Gran+", 
                         paste0("PC",1:20,collapse="+"))
covars_mat <- as.matrix(model.frame(as.formula(covars_formula), nonMethData, na.action=na.pass))
complete_cases <- complete.cases(covars_mat)
design_mat <- cbind(covars_mat, t(Bvals))[complete_cases,]

# Additional needed variables
shareids <- nonMethData$shareid[complete_cases]
cvd_surv <- Surv(time=nonMethData$timeToEvent, event=nonMethData$event)[complete_cases]  # Outcome
frs <- calc_FRS(nonMethData)[complete_cases]

trainTestMRS <- function(cvd_surv, design_mat, frs, shareids, trainset_idx) {
  train_ids <- shareids[trainset_idx]
  mrs.fit <- glmnet(design_mat[trainset_idx,], cvd_surv[trainset_idx], family="cox", alpha=0.5,
                    penalty.factor=ifelse((grepl("^cg", colnames(design_mat)) | 
                                             grepl("^PC", colnames(design_mat))), 1, 0))
  coefs <- mrs.fit$beta[,ncol(mrs.fit$beta)]
  coefs_meth <- coefs[coefs!=0 & grepl("cg", names(coefs))]
  mrs_calc <- as.vector(design_mat[,names(coefs_meth)] %*% coefs_meth)
  testset <- setdiff(1:nrow(design_mat), trainset_idx)
  mrs.test <- coxph(cvd_surv~mrs+frs, data=data.frame(cvd_surv=cvd_surv, mrs=mrs_calc, frs=frs), subset=testset)
  list(components=coefs_meth, res=summary(mrs.test)$coef, train_ids=train_ids)
}

allRes <- list()
set.seed(1)
for (i in 1:20) {
  print(paste("Iteration", i))
  twoFold <- createFolds(factor(cvd_surv[,"status"]), 2)
  allRes[[i]] <- trainTestMRS(cvd_surv, design_mat, frs, shareids, twoFold[[1]])
}

save("allRes", file="../int/testMRS_repCV_res_BVALS_lowPen.RData")


