suppressMessages(silent <- lapply(c("tidyverse","doParallel","itertools","survival","gtools","coxme","glmnet",
                                    "caret","Coxnet"), library, character.only=T))

# Load methylation and sample data
load("../int/Mvals.RData")  # Loads Mvals object
print("Methylation data loaded.")

# Load non-methylation (event + covariate) data and sample similarity matrix
load("../int/nonMethData.RData")
nonMethData$survObj <- Surv(time=nonMethData$timeToEvent, event=nonMethData$event, type="right")
load("../int/Kmat.RData")

# Sanity checks
final_id_set <- intersect(nonMethData$sampleKey, colnames(Mvals))
nonMethData <- nonMethData[match(final_id_set, nonMethData$sampleKey),]
Mvals <- Mvals[,final_id_set]
Kmat <- Kmat[final_id_set, final_id_set]
print(paste("Dimensions of event/covariate matrix:", dim(nonMethData)))
print(paste("Dimensions of M-value matrix:", dim(Mvals)))
stopifnot(ncol(Mvals)==nrow(nonMethData),  # Ensure same # of samples in methylation and covariate data
          all(colnames(Mvals)==nonMethData$sampleKey),  # Ensure identical order of samples for methylation and covariate data
          all(colnames(Mvals)==rownames(Kmat)))  # Ensure similarity matrix is in the correct order  

## Run regressions

traintest_MRS <- function(test_idx, methData, nonMethData) {
  
  myTry <- function(expr, CpG) {
    # Captures model failures and returns successful result or vector of NAs
    tryCatch(expr,
             error = function(e) {print(e); return(c(CpG, rep(NA, 2)))},
             warning = function(w) {print(w); return(c(CpG, rep(NA, 2)))})
  }
  
  # Formulas for fixed-effects only models of interest
  run_cox <- function(probeData, model_spec) {
    # Given a row of the M-value matrix (corresponding to a CpG site), bind that methylation
    # data to the covariate data and run Cox proportional hazards regression
    CpG <- rownames(probeData)
    modelData <- cbind(trainNonMethData, meth=as.numeric(probeData))
    myTry({
      cox.fit <- coxph(as.formula(model_spec), data=modelData)
      c(CpG=CpG, summary(cox.fit)$coef['meth',c('coef','Pr(>|z|)')])
      # lr.fit <- glm(as.formula(model_spec), family="binomial", data=modelData)
      # c(CpG=CpG, summary(lr.fit)$coefficients['meth',c('Estimate','Pr(>|z|)')])
    }, CpG)
  }
  
  train_MRS <- function(train_idx, modelData, penalized=F) {
    
    trainSet <- modelData[train_idx,]
    
    if (penalized==T) {
      design_mat_form <- paste0("~sex+age+smoking_now+bmi+CD8T+CD4T+NK+Bcell+Mono+Gran+", 
                                paste0("PC",1:20,"_cp",collapse="+"),"+",paste0(keep_CpGs,collapse="+"))
      x <- model.frame(as.formula(design_mat_form), trainSet, na.action=na.pass)
      x <- apply(x, 2, scale)
      y <- cbind(time=trainSet$timeToEvent, status=trainSet$event)
      # y <- trainSet$event
      mrs.fit <- Coxnet(x, y, penalty="Enet", alpha=0.5)
      # mrs.fit <- glmnet(x, y, family="binomial", alpha=0.5)
      coefs <- mrs.fit$Beta[which(colnames(x) %in% keep_CpGs), ncol(mrs.fit$Beta)]
      # coefs <- mrs.fit$beta[which(colnames(x) %in% keep_CpGs), ncol(mrs.fit$beta)]
    } else {
      model_spec <- paste0("survObj~sex+age+smoking_now+bmi+CD8T+CD4T+NK+Bcell+Mono+Gran+", 
                           paste0("PC",1:20,"_cp",collapse="+"),"+",paste0(keep_CpGs,collapse="+"))
      mrs.fit <- coxph(as.formula(model_spec), data=trainSet)
      coefs <- summary(mrs.fit)$coef[keep_CpGs,'coef']  # Extract regression coefficients
    }
    coefs
  }
  
  test_MRS <- function(test_idx, modelData, coefs, keep_CpGs) {
    testSet <- modelData[test_idx,]
    testSet$mrs <- as.vector(as.matrix(testSet[,keep_CpGs]) %*% coefs)
    
    source("helpers.R")
    testSet$frs <- calc_FRS(testSet[,-which(colnames(testSet)=="survObj")])
    
    FRSonly.fit <- coxph(survObj~frs, data=testSet)
    FRSonly.res <- summary(FRSonly.fit)$coef['frs',c(1,5)]
    MRSonly.fit <- coxph(survObj~mrs, data=testSet)
    MRSonly.res <- summary(MRSonly.fit)$coef['mrs',c(1,5)]
    combined.fit <- coxph(survObj~frs+mrs, data=testSet)
    # combined.fit <- glm(event~frs+mrs, data=testSet)
    combined.res <- summary(combined.fit)$coef[c('frs','mrs'),c(1,5)]
    # combined.res <- summary(combined.fit)$coefficients[c('frs','mrs'),c(1,5)]
    res <- c(FRSonly.res, MRSonly.res, combined.res['frs',], combined.res['mrs',])
    # res <- combined.res
    names(res) <- paste0(rep(c("FRS","MRS"), each=2, length.out=8),
                         rep(c("coef","p"), length.out=8),
                         rep(c("alone","adjusted"), each=4))
    res
  }
  
  # run_MRS <- function(test_idx, modelData, penalized=F) {
  #   ## Creates an MRS based on an initial regression, then tests compared to the FRS
  #   coefs <- train_MRS(setDiff(1:nrow(modelData), test_idx), penalized)
  #   test_MRS(test_idx, modelData, coefs)
  # }
  
  model_list <- list(
    basic_wbc_20PC=paste0("survObj~meth+sex+age+smoking_now+bmi+CD8T+CD4T+NK+Bcell+Mono+Gran+", 
                          paste0("PC",1:20,collapse="+"))
  )
  model_list_lr <- list(
    basic_wbc_20PC=paste0("event~meth+sex+age+smoking_now+bmi+CD8T+CD4T+NK+Bcell+Mono+Gran+", 
                          paste0("PC",1:20,collapse="+"))
  )
  
  # Determine train and test sets
  train_idx <- setdiff(1:nrow(nonMethData), test_idx)
  trainMvals <- methData[,train_idx]
  trainNonMethData <- nonMethData[train_idx,]
  # Parallelized Cox EWAS for incident CVD
  cl <- makePSOCKcluster(detectCores(), outfile="progress.txt")
  registerDoParallel(cl)
  res <- foreach(meth=enumerate(iter(trainMvals, by="row")), .combine=rbind, .packages=c("survival")) %dopar%
    run_cox(meth$value, model_list[["basic_wbc_20PC"]])
  stopImplicitCluster()
  # Determine set of candidate CpGs
  resDF <- data.frame(res, stringsAsFactors=F)
  names(resDF) <- c("CpG","beta","p")
  resDF$p <- as.numeric(resDF$p)
  resDF$fdr <- p.adjust(resDF$p, method="BH")
  keep_CpGs <- na.omit(resDF$CpG[resDF$fdr < 0.05])
  print(head(keep_CpGs))
  # Regression to create MRS
  Mvals_keep <- Mvals[keep_CpGs,]
  allModelData <- cbind(nonMethData, t(Mvals_keep))
  coefs <- train_MRS(train_idx, allModelData, penalized=T)
  testRes <- test_MRS(test_idx, allModelData, coefs, keep_CpGs)
  print(testRes)
  testRes
}

set.seed(1)
folds <- createFolds(nonMethData$event, 3)  # For 3-fold CV, keep proportions of people with events constant
lapply(folds, traintest_MRS, Mvals, nonMethData)
# train_idx <- setdiff(1:nrow(nonMethData), folds[[1]])
# test_idx <- folds[[1]]
# trainMvals <- Mvals[,train_idx]
# trainNonMethData <- nonMethData[train_idx,]
# res <- foreach(meth=enumerate(iter(trainMvals, by="row")), .combine=rbind, .packages=c("survival")) %dopar% {
#     run_cox(meth$value, model_list[["basic_wbc_20PC"]])
#   }

# 
# ## Determine set of candidate CpGs from EWAS results
# resDF <- data.frame(res, stringsAsFactors=F)
# names(resDF) <- c("CpG","beta","p")
# resDF$p <- as.numeric(resDF$p)
# resDF$fdr <- p.adjust(resDF$p, method="BH")
# keep_CpGs <- na.omit(resDF$CpG[resDF$fdr < 0.05])
# print(head(keep_CpGs))
# 
# Mvals_keep <- Mvals[keep_CpGs,final_id_set]
# allModelData <- cbind(nonMethData, t(Mvals_keep))
# 
# 
# coefs <- train_MRS(train_idx, allModelData, penalized=T)
# test <- allModelData[test_idx,]
# print(test_MRS(test_idx, test, coefs))

# # Generate K-M plot
# test$mrs <- as.vector(as.matrix(test[,keep_CpGs]) %*% coefs)
# test$mrs_binary <- test$mrs > mean(test$mrs)
# test.survfit <- survfit(survObj~mrs_binary, data=test, conf.type="log-log")
# test.survfit
# jpeg("../output/survplot.jpg")
# plot(test.survfit, col=1:2)
# legend("left", legend=c("lowMRS","highMRS"), col=1:2, lty=1)
# dev.off()

