suppressMessages(silent <- lapply(c("tidyverse","doParallel","itertools","survival","gtools","coxme","glmnet","minfi","caret"), 
                                  library, character.only=T))
source("helpers.R")

### CREATE MRS

# Load methylation and sample data
load("../int/Mvals.RData")  # Loads Mvals object
print("Methylation data loaded.")

# Load non-methylation (event + covariate) data and sample similarity matrix
load("../int/nonMethData2.RData")

# Sanity checks
final_id_set <- intersect(nonMethData$sampleKey, colnames(Mvals))
nonMethData <- nonMethData[match(final_id_set, nonMethData$sampleKey),]
Mvals <- Mvals[,final_id_set]
print(paste("Dimensions of event/covariate matrix:", dim(nonMethData)))
print(paste("Dimensions of M-value matrix:", dim(Mvals)))
stopifnot(ncol(Mvals)==nrow(nonMethData),  # Ensure same # of samples in methylation and covariate data
          all(colnames(Mvals)==nonMethData$sampleKey))  # Ensure identical order of samples for methylation and covariate data

# Covariates for regressions (including methylation M-values)
covars_formula <- paste0("~sex+age+smoking_now+bmi+CD8T+CD4T+NK+Bcell+Mono+Gran+", 
                          paste0("PC",1:20,collapse="+"))
covars_mat <- as.matrix(model.frame(as.formula(covars_formula), nonMethData, na.action=na.pass))
design_mat <- cbind(covars_mat, t(Mvals))

# Dependent variables for regressions
cvd <- Surv(time=nonMethData$timeToEvent, event=nonMethData$event)
frs <- calc_FRS(nonMethData)
imt <- nonMethData$IMT_mean

trainTestMRS <- function(dep_var, design_mat, cvd, frs, testset, type) {
  # complete <- complete.cases(dep_var) & complete.cases(design_mat)
  print(type)
  trainset <- setdiff(1:nrow(design_mat), testset)
  mrs.fit <- glmnet(design_mat[trainset,], dep_var[trainset], family=type, alpha=0.5)
  coefs <- mrs.fit$beta[,ncol(mrs.fit$beta)]
  coefs_meth <- coefs[coefs!=0 & grepl("cg", names(coefs))]
  mrs_calc <- as.vector(design_mat[,names(coefs_meth)] %*% coefs_meth)
  mrs.test <- coxph(cvd~mrs+frs, data=data.frame(cvd=cvd, mrs=mrs_calc, frs=frs), subset=testset)
  # print(summary(mrs.test)$coef)
  list(num_cpgs=sum(coefs!=0), res=summary(mrs.test)$coef)
}

# complete_subjects <- complete.cases(design_mat) & !is.na(cvd) & !is.na(frs) & !is.na(imt)
outcomes <- list(cvd=list(dep_var=cvd, type="cox"),
     frs=list(dep_var=frs, type="gaussian"),
     imt=list(dep_var=imt, type="gaussian"))
lapply(outcomes, function(outcome) {
  set.seed(100)
  complete_subjects <- !is.na(outcome$dep_var) & complete.cases(design_mat)
  folds <- createFolds(factor(cvd[complete_subjects,"status"]), 3)
  lapply(folds, function(testset) {
    trainTestMRS(outcome$dep_var[complete_subjects], design_mat[complete_subjects,], 
                 cvd[complete_subjects], frs[complete_subjects],
                 testset, outcome$type)
  })
})

