suppressMessages(silent <- lapply(c("tidyverse","doParallel","itertools","survival"), 
                                  library, character.only=T))

# Load methylation data
load("../int/Mvals.RData")  # Loads Mvals object
print("Methylation data loaded.")

# Load non-methylation (event + covariate) data
load("../int/nonMethData.RData")
nonMethData$survObj <- Surv(time=nonMethData$timeToEvent, event=nonMethData$event, type="right")

# Sanity checks
Mvals <- Mvals[,match(nonMethData$sampleKey, colnames(Mvals))]
print(paste("Dimensions of event/covariate matrix:", paste(dim(nonMethData), collapse=" x ")))
print(paste("Dimensions of M-value matrix:", paste(dim(Mvals), collapse=" x ")))
stopifnot(all(colnames(Mvals)==nonMethData$sampleKey))  # Ensure identical number and order of samples for methylation and covariate data

## Run regressions
set.seed(1)
Mvals <- Mvals[sample(1:nrow(Mvals), 100000),]
cl <- makePSOCKcluster(detectCores(), outfile="progress.txt")
registerDoParallel(cl)

myTry <- function(expr, CpG) {
  # Captures model failures and returns successful result or vector of NAs
  tryCatch(expr,
           error = function(e) {print(e); return(c(CpG, rep(NA, 3)))},
           warning = function(w) {print(w); return(c(CpG, rep(NA, 3)))})
}

# Formulas for fixed-effects only models of interest
run_cox <- function(probeData, model_spec) {
  # Given a row of the M-value matrix (corresponding to a CpG site), bind that methylation
  # data to the covariate data and run Cox proportional hazards regression
  CpG <- rownames(probeData)
  modelData <- cbind(nonMethData, meth=as.numeric(probeData))
  myTry({
    cox.fit <- coxph(as.formula(model_spec), data=modelData)
    c(CpG=CpG, summary(cox.fit)$coef['meth',c('coef','z','Pr(>|z|)')])
  }, CpG)
} 

model_list <- list(
  basic=paste0("survObj~meth+sex+age+smoking_now+bmi"),
  wbcOnly=paste0("survObj~meth+sex+age+smoking_now+bmi+CD4T+NK+Bcell+Mono+Gran"),
  pastEventOnly=paste0("survObj~meth+sex+age+smoking_now+bmi+pastEvent"),
  wbcPastEvent=paste0("survObj~meth+sex+age+smoking_now+bmi+CD4T+NK+Bcell+Mono+Gran+pastEvent"),
  wbcPastEvent1CPA=paste0("survObj~meth+sex+age+smoking_now+bmi+CD4T+NK+Bcell+Mono+Gran+pastEvent+PC1_cp")
)

res <- foreach(mod=model_list, .final=function(l) setNames(l, names(model_list))) %do% {
  foreach(meth=enumerate(iter(Mvals, by="row")), .combine=rbind, .packages=c("survival")) %dopar% {
    if(meth$index %% 10000 == 0) print(paste0("Model ", mod, ", CpG ", meth$index))
    run_cox(meth$value, mod)
  }
}

stopCluster(cl)

print("Saving results...")
save("res", file=paste0("../int/ewasRes7_pastEventsAdjustment.RData"))