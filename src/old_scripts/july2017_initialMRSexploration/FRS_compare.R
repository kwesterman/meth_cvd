library(tidyverse)
library(survival)
library(gtools)
library(caret)
library(Coxnet)

## Determine set of candidate CpGs from EWAS results
load("../int/ewasRes.RData")
resDF <- data.frame(res, stringsAsFactors=F)
names(resDF) <- c("CpG","beta","p")
resDF$p <- as.numeric(resDF$p)
resDF$fdr <- p.adjust(resDF$p, method="BH")
# keep_CpG_idx <- order(resDF$fdr)[1:20]
keep_CpGs <- na.omit(resDF$CpG[resDF$fdr < 0.05])
keep_CpGs <- resDF$CpG[keep_CpG_idx]
print(head(keep_CpGs))

# Load M-values
load("../int/Mvals.RData")

# Load non-methylation (event + covariate) data
load("../int/nonMethData.RData")
nonMethData$survObj <- Surv(time=nonMethData$timeToEvent, event=nonMethData$event, type="right")

final_id_set <- intersect(nonMethData$sampleKey, colnames(Mvals))
nonMethData <- nonMethData[match(final_id_set, nonMethData$sampleKey),]
Mvals_keep <- Mvals[keep_CpGs,final_id_set]
allModelData <- cbind(nonMethData, t(Mvals_keep))

train_MRS <- function(train_idx, modelData, penalized=F) {
  
  trainSet <- allModelData[train_idx,]
  if (penalized==T) {
    design_mat_form <- paste0("~sex+age+smoking_now+bmi+CD8T+CD4T+NK+Bcell+Mono+Gran+", 
                              paste0("PC",1:20,"_cp",collapse="+"),"+",paste0(keep_CpGs,collapse="+"))
    x <- model.frame(as.formula(design_mat_form), trainSet, na.action=na.pass)
    x <- apply(x, 2, scale)
    y <- cbind(time=trainSet$timeToEvent, status=trainSet$event)
    mrs.fit <- Coxnet(x, y, penalty="Enet", alpha=0.5)
    testSet <- allModelData[test_idx,]
    coefs <- mrs.fit$Beta[which(colnames(x) %in% keep_CpGs), ncol(mrs.fit$Beta)]
  } else {
    model_spec <- paste0("survObj~sex+age+smoking_now+bmi+CD8T+CD4T+NK+Bcell+Mono+Gran+", 
                         paste0("PC",1:20,"_cp",collapse="+"),"+",paste0(keep_CpGs,collapse="+"))
    mrs.fit <- coxph(as.formula(model_spec), data=trainSet)
    coefs <- summary(mrs.fit)$coef[keep_CpGs,'coef']  # Extract regression coefficients
  }
  coefs
}

test_MRS <- function(test_idx, modelData, coefs) {
  testSet <- allModelData[test_idx,]
  testSet$mrs <- as.vector(as.matrix(testSet[,keep_CpGs]) %*% coefs)
  
  source("helpers.R")
  testSet$frs <- calc_FRS(testSet[,-which(colnames(testSet)=="survObj")])
  
  FRSonly.fit <- coxph(survObj~frs, data=testSet)
  FRSonly.res <- summary(FRSonly.fit)$coef['frs',c(1,5)]
  MRSonly.fit <- coxph(survObj~mrs, data=testSet)
  MRSonly.res <- summary(MRSonly.fit)$coef['mrs',c(1,5)]
  combined.fit <- coxph(survObj~frs+mrs, data=testSet)
  combined.res <- summary(combined.fit)$coef[c('frs','mrs'),c(1,5)]
  res <- c(FRSonly.res, MRSonly.res, combined.res['frs',], combined.res['mrs',])
  names(res) <- paste0(rep(c("FRS","MRS"), each=2, length.out=8), 
                       rep(c("coef","p"), length.out=8),
                       rep(c("alone","adjusted"), each=4))
  res
}

run_MRS <- function(test_idx, modelData, penalized=F) {
  ## Creates an MRS based on an initial regression, then tests compared to the FRS
  coefs <- train_MRS(setdiff(1:nrow(modelData), test_idx), modelData, penalized)
  test_MRS(test_idx, modelData, coefs)
}

set.seed(1)
folds <- createFolds(allModelData$event, 5)  # For 5-fold CV, keeping proportions of people with events constant
modelStats <- do.call(rbind, lapply(folds, run_MRS, allModelData, penalized=F))
summary(modelStats)
modelStats_penalized <- do.call(rbind, lapply(folds, run_MRS, allModelData, penalized=T))
summary(modelStats_penalized)

folds <- createFolds(allModelData$event, 3)
coefs <- train_MRS(setdiff(1:nrow(allModelData), folds[[1]]), allModelData, penalized=T)
test <- allModelData[test_idx,]
test$mrs <- as.vector(as.matrix(test[,keep_CpGs]) %*% coefs)
test$mrs_binary <- test$mrs > mean(test$mrs)
test.survfit <- survfit(survObj~mrs_binary, data=test, conf.type="log-log")
test.survfit
plot(test.survfit, col=1:2)
legend("left", legend=c("lowMRS","highMRS"), col=1:2, lty=1)
