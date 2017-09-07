suppressMessages(silent <- lapply(c("tidyverse","doParallel","itertools"), 
                                  library, character.only=T))

numCores <- 100

# Load methylation data
load("../int/Mvals.RData")  # Loads Mvals object
print("Methylation data loaded.")

# Load non-methylation (event + covariate) data
load("../int/nonMethData.RData")

# Sanity checks
Mvals <- Mvals[,match(nonMethData$sampleKey, colnames(Mvals))]
print(paste("Dimensions of event/covariate matrix:", paste(dim(nonMethData), collapse=" x ")))
print(paste("Dimensions of M-value matrix:", paste(dim(Mvals), collapse=" x ")))
stopifnot(all(colnames(Mvals)==nonMethData$sampleKey))  # Ensure identical number and order of samples for methylation and covariate data

## Run regressions
# Mvals <- Mvals[sample(1:nrow(Mvals), 50000),]
cl <- makePSOCKcluster(numCores, outfile="progress.txt")
registerDoParallel(cl)

myTry <- function(expr, CpG) {
  # Captures model failures and returns successful result or vector of NAs
  tryCatch(expr,
           error = function(e) {print(e); return(c(CpG, rep(NA, 3)))},
           warning = function(w) {print(w); return(c(CpG, rep(NA, 3)))})
}

# Formulas for fixed-effects only models of interest
run_logreg <- function(probeData, model_spec) {
  # Given a row of the M-value matrix (corresponding to a CpG site), bind that methylation
  # data to the covariate data and run logreg proportional hazards regression
  CpG <- rownames(probeData)
  modelData <- cbind(nonMethData, meth=as.numeric(probeData))
  myTry({
    logreg.fit <- glm(as.formula(model_spec), data=modelData, family="binomial")
    c(CpG=CpG, summary(logreg.fit)$coef['meth',c('Estimate','z value','Pr(>|z|)')])
  }, CpG)
} 

model_list <- list(
  sixPCs=paste0("event~meth+sex+age+smoking_now+bmi+", paste0("PC", 1:6, collapse="+")),
  fourPCsWbc=paste0("event~meth+sex+age+smoking_now+bmi+CD8T+CD4T+NK+Bcell+Mono+Gran+", paste0("PC", 1:4, collapse="+"))
)

res <- foreach(mod=model_list, .final=function(l) setNames(l, names(model_list))) %do% {
  foreach(meth=enumerate(iter(Mvals, by="row")), .combine=rbind) %dopar% {
    if(meth$index %% 10000 == 0) print(paste0("Model ", mod, ", CpG ", meth$index))
    run_logreg(meth$value, mod)
  }
}

stopCluster(cl)

print("Saving results...")
save("res", file=paste0("../int/lrEwasRes.RData"))